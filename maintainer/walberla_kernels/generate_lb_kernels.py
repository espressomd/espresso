#
# Copyright (C) 2020-2021 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import os
import sys
import jinja2
import argparse

import sympy as sp
import pystencils as ps
from lbmpy.creationfunctions import create_lb_collision_rule, create_mrt_orthogonal
import lbmpy.forcemodels
from pystencils_walberla import CodeGeneration, codegen

from lbmpy_walberla import generate_boundary, generate_lb_pack_info
from walberla_lbm_generation import generate_macroscopic_values_accessors

import lbmpy.boundaries
import lbmpy_walberla.additional_data_handler
import relaxation_rates
from lbmpy.fieldaccess import CollideOnlyInplaceAccessor
from lbmpy.stencils import get_stencil
from lbmpy.updatekernels import create_lbm_kernel, create_stream_pull_with_output_kernel
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter

import lees_edwards
# for collide-push from lbmpy.fieldaccess import StreamPushTwoFieldsAccessor

from lbmpy.advanced_streaming.indexing import NeighbourOffsetArrays


parser = argparse.ArgumentParser(description='Generate the waLBerla kernels.')
parser.add_argument('--single-precision', action='store_true', required=False,
                    help='Use single-precision')
args = parser.parse_args()

data_type_cpp = {True: 'double', False: 'float'}
data_type_np = {True: 'float64', False: 'float32'}


def adapt_pystencils():
    '''
    Adapt pystencils to the SFINAE method (add the block offset lambda
    callback and the time_step increment).
    '''
    old_add_pystencils_filters_to_jinja_env = codegen.add_pystencils_filters_to_jinja_env

    def new_add_pystencils_filters_to_jinja_env(jinja_env):
        # save original pystencils to adapt
        old_add_pystencils_filters_to_jinja_env(jinja_env)
        old_generate_members = jinja_env.filters['generate_members']
        old_generate_refs_for_kernel_parameters = jinja_env.filters[
            'generate_refs_for_kernel_parameters']

        @jinja2.contextfilter
        def new_generate_members(*args, **kwargs):
            output = old_generate_members(*args, **kwargs)
            token = ' block_offset_0_;'
            if token in output:
                i = output.index(token)
                vartype = output[:i].split('\n')[-1].strip()
                output += f'\nstd::function<void(IBlock *, {vartype}&, {vartype}&, {vartype}&)> block_offset_generator = [](IBlock * const, {vartype}&, {vartype}&, {vartype}&) {{ }};'
            return output

        def new_generate_refs_for_kernel_parameters(*args, **kwargs):
            output = old_generate_refs_for_kernel_parameters(*args, **kwargs)
            if 'block_offset_0' in output:
                old_token = 'auto & block_offset_'
                new_token = 'auto block_offset_'
                assert output.count(old_token) == 3, \
                    f'could not find "{old_token}" in """\n{output}\n"""'
                output = output.replace(old_token, new_token)
                output += '\nblock_offset_generator(block, block_offset_0, block_offset_1, block_offset_2);'
            return output

        # replace pystencils
        jinja_env.filters['generate_members'] = new_generate_members
        jinja_env.filters['generate_refs_for_kernel_parameters'] = new_generate_refs_for_kernel_parameters

    codegen.add_pystencils_filters_to_jinja_env = new_add_pystencils_filters_to_jinja_env


def earmark_generated_kernels():
    '''
    Add an earmark at the beginning of generated kernels to document the
    pystencils/lbmpy toolchain that was used to create them.
    '''
    walberla_root = lbmpy_walberla.__file__.split('/python/lbmpy_walberla/')[0]
    with open(os.path.join(walberla_root, '.git/HEAD')) as f:
        walberla_commit = f.read()
    if walberla_commit.startswith('ref: refs/heads/master'):
        ref = walberla_commit.split()[1]
        with open(os.path.join(walberla_root, f'.git/{ref}')) as f:
            walberla_commit = f.read()
    token = '// kernel generated with'
    earmark = (
        f'{token} pystencils v{ps.__version__}, lbmpy v{lbmpy.__version__}, '
        f'lbmpy_walberla/pystencils_walberla from commit {walberla_commit}\n'
    )
    for filename in os.listdir('.'):
        if filename.endswith(('.h', '.cpp')):
            with open(filename, 'r+') as f:
                content = f.read()
                if not content.startswith(token):
                    f.seek(0)
                    f.write(earmark + content)


def guard_generated_kernels_clang_format():
    '''
    Some namespaces are too long and will break ``clang-format`` versions
    9 and 10. Replace them with a unique string of reasonable size.
    '''
    import re
    import hashlib
    for filename in os.listdir('.'):
        if filename.endswith('.cpp'):
            with open(filename, 'r') as f:
                content = f.read()
            all_ns = re.findall(r"^namespace (internal_[a-zA-Z0-9_]{54,}) \{$",
                                content, flags=re.MULTILINE)
            if not all_ns:
                continue
            for ns in all_ns:
                content = re.sub(rf"(?<=[^a-zA-Z0-9_]){ns}(?=[^a-zA-Z0-9_])",
                                 f"internal_{hashlib.md5(ns.encode('utf-8')).hexdigest()}",
                                 content)
            with open(filename, 'w') as f:
                f.write(content)


def generate_fields(ctx, stencil):
    dtype = data_type_np[ctx.double_accuracy]
    field_layout = 'fzyx'
    q = len(stencil)
    dim = len(stencil[0])

    fields = {}
    # Symbols for PDF (twice, due to double buffering)
    fields['pdfs'] = ps.Field.create_generic(
        'pdfs',
        dim,
        dtype,
        index_dimensions=1,
        layout=field_layout,
        index_shape=(q,)
    )
    fields['pdfs_tmp'] = ps.Field.create_generic(
        'pdfs_tmp',
        dim,
        dtype,
        index_dimensions=1,
        layout=field_layout,
        index_shape=(q,)
    )
    fields['velocity'] = ps.Field.create_generic(
        'velocity',
        dim,
        dtype,
        index_dimensions=1,
        layout=field_layout,
        index_shape=(dim,)
    )
    fields['force'] = ps.Field.create_generic(
        'force',
        dim,
        dtype,
        index_dimensions=1,
        layout=field_layout,
        index_shape=(dim,)
    )

    return fields


def generate_collision_sweep(
        ctx, lb_method, collision_rule, class_name, params):

    # Symbols for PDF (twice, due to double buffering)
    fields = generate_fields(ctx, lb_method.stencil)

    # Generate collision kernel
    collide_update_rule = create_lbm_kernel(
        collision_rule,
        fields['pdfs'],
        fields['pdfs_tmp'],
        CollideOnlyInplaceAccessor())
    collide_ast = ps.create_kernel(collide_update_rule, **params,
                                   data_type=data_type_np[ctx.double_accuracy])
    collide_ast.function_name = 'kernel_collide'
    collide_ast.assumed_inner_stride_one = True
    codegen.generate_sweep(ctx, class_name, collide_ast)


def generate_stream_sweep(ctx, lb_method, class_name, params):
    # Symbols for PDF (twice, due to double buffering)
    fields = generate_fields(ctx, lb_method.stencil)

    # Generate stream kernel
    stream_update_rule = create_stream_pull_with_output_kernel(lb_method, fields['pdfs'], fields['pdfs_tmp'],
                                                               output={'velocity': fields['velocity']})
    stream_ast = ps.create_kernel(stream_update_rule, **params,
                                  data_type=data_type_np[ctx.double_accuracy])
    stream_ast.function_name = 'kernel_stream'
    stream_ast.assumed_inner_stride_one = True
    codegen.generate_sweep(
        ctx, class_name, stream_ast, field_swaps=[(fields['pdfs'], fields['pdfs_tmp'])])


def generate_setters(lb_method):
    fields = generate_fields(ctx, lb_method.stencil)

    initial_rho = sp.Symbol('rho_0')
    pdfs_setter = macroscopic_values_setter(lb_method,
                                            initial_rho,
                                            fields['velocity'].center_vector,
                                            fields['pdfs'].center_vector)
    return pdfs_setter


def check_dependencies():
    import setuptools
    import pystencils
    import lbmpy
    SpecifierSet = setuptools.version.pkg_resources.packaging.specifiers.SpecifierSet
    for module, requirement in [(pystencils, '==0.4.4'), (lbmpy, '==0.4.4')]:
        assert SpecifierSet(requirement).contains(module.__version__), \
            f"{module.__name__} version {module.__version__} doesn't match requirement {requirement}"


class BounceBackSlipVelocityUBB(
        lbmpy_walberla.additional_data_handler.UBBAdditionalDataHandler):
    '''
    Dynamic UBB that implements the bounce-back method with slip velocity.
    '''

    def data_initialisation(self, direction):
        '''
        Modified ``indexVector`` initialiser. The "classical" dynamic UBB
        uses the velocity callback as a velocity flow profile generator.
        Here we use that callback as a bounce-back slip velocity generator.
        This way, the dynamic UBB can be used to implement a LB boundary.
        '''
        code = super().data_initialisation(direction)
        dirVec = self.stencil_info[direction][1]
        token = ' = elementInitaliser(Cell(it.x(){}, it.y(){}, it.z(){}),'
        old_initialiser = token.format('', '', '')
        assert old_initialiser in code
        new_initialiser = token.format(
            '+' + str(dirVec[0]),
            '+' + str(dirVec[1]),
            '+' + str(dirVec[2])).replace('+-', '-')
        return code.replace(old_initialiser, new_initialiser)


class PatchedUBB(lbmpy.boundaries.UBB):
    '''
    Velocity bounce back boundary condition, enforcing specified velocity at obstacle.
    '''

    def __call__(self, f_out, f_in, dir_symbol,
                 inv_dir, lb_method, index_field):
        '''
        Modify the assignments such that the source and target pdfs are swapped.
        '''
        assignments = super().__call__(
            f_out, f_in, dir_symbol, inv_dir, lb_method, index_field)

        assert len(assignments) > 0

        out = []
        if len(assignments) > 1:
            out.extend(assignments[:-1])

        neighbor_offset = NeighbourOffsetArrays.neighbour_offset(
            dir_symbol, lb_method.stencil)

        assignment = assignments[-1]
        assert assignment.lhs.field == f_in
        out.append(ps.Assignment(assignment.lhs.get_shifted(*neighbor_offset),
                                 assignment.rhs - f_out(dir_symbol) + f_in(dir_symbol)))
        return out


check_dependencies()
adapt_pystencils()


class PatchedCodeGeneration(CodeGeneration):

    def __init__(self):
        old_sys_argv = sys.argv
        sys.argv = sys.argv[:1]
        super().__init__()
        sys.argv = old_sys_argv


with PatchedCodeGeneration() as ctx:
    ctx.double_accuracy = not args.single_precision
    precision_prefix = {
        True: 'DoublePrecision',
        False: 'SinglePrecision'}[
        ctx.double_accuracy]
    precision_suffix = {
        True: 'double_precision',
        False: 'single_precision'}[
        ctx.double_accuracy]
    precision_rng = {
        True: ps.rng.PhiloxTwoDoubles,
        False: ps.rng.PhiloxFourFloats}[
        ctx.double_accuracy]
    kT = sp.symbols('kT')
    stencil = get_stencil('D3Q19')
    fields = generate_fields(ctx, stencil)
    force_field = fields['force']

    # Vectorization parameters
    cpu_vectorize_info = {
        "instruction_set": "avx",
        "assume_inner_stride_one": True,
        "assume_aligned": True,
        "assume_sufficient_line_padding": False}
    params = {"target": "cpu"}
    params_vec = {"target": "cpu", "cpu_vectorize_info": cpu_vectorize_info}

    # LB Method definition
    method = create_mrt_orthogonal(
        stencil=stencil,
        compressible=True,
        weighted=True,
        relaxation_rates=relaxation_rates.rr_getter,
        force_model=lbmpy.forcemodels.Schiller(force_field.center_vector)
    )

    # generate stream kernels
    generate_stream_sweep(
        ctx,
        method,
        f"StreamSweep{precision_prefix}",
        params)
    generate_stream_sweep(
        ctx,
        method,
        f"StreamSweep{precision_prefix}AVX",
        params_vec)

    # generate initial densities
    pdfs_setter = generate_setters(method)
    codegen.generate_sweep(
        ctx,
        f"InitialPDFsSetter{precision_prefix}",
        pdfs_setter)

    # generate unthermalized collision rule
    collision_rule_unthermalized = create_lb_collision_rule(
        method,
        optimization={'cse_global': True,
                      'double_precision': ctx.double_accuracy}
    )
    generate_collision_sweep(
        ctx,
        method,
        collision_rule_unthermalized,
        f"CollideSweep{precision_prefix}",
        {}
    )
    generate_collision_sweep(
        ctx,
        method,
        collision_rule_unthermalized,
        f"CollideSweep{precision_prefix}AVX",
        {"cpu_vectorize_info": cpu_vectorize_info}
    )

    # generate unthermalized Lees-Edwards collision rule

    le_config = lbmpy.LBMConfig(stencil=stencil, method=lbmpy.Method.TRT, relaxation_rate=sp.Symbol("omega_shear"), compressible=True,
                                force_model=lbmpy.ForceModel.GUO,
                                force=force_field.center_vector, kernel_type='collide_only')
    lbm_opt = lbmpy.LBMOptimisation(symbolic_field=fields["pdfs"])
    le_collision_rule_unthermalized = lbmpy.create_lb_update_rule(
        lbm_config=le_config,
        lbm_optimisation=lbm_opt)

    le_collision_rule_unthermalized = lees_edwards.add_lees_edwards_to_collision(
        le_collision_rule_unthermalized,
        fields["pdfs"], stencil, 1)  # shear_dir_normal y
    generate_collision_sweep(
        ctx,
        le_config,
        le_collision_rule_unthermalized,
        f"CollideSweep{precision_prefix}LeesEdwards",
        {}
    )
    generate_collision_sweep(
        ctx,
        le_config,
        le_collision_rule_unthermalized,
        f"CollideSweep{precision_prefix}LeesEdwardsAVX",
        {"cpu_vectorize_info": cpu_vectorize_info}
    )

    # generate thermalized LB
    collision_rule_thermalized = create_lb_collision_rule(
        method,
        fluctuating={
            'temperature': kT,
            'block_offsets': 'walberla',
            'rng_node': precision_rng
        },
        optimization={'cse_global': True,
                      'double_precision': ctx.double_accuracy}
    )
    generate_collision_sweep(
        ctx,
        method,
        collision_rule_thermalized,
        f"CollideSweep{precision_prefix}Thermalized",
        params
    )
    generate_collision_sweep(
        ctx,
        method,
        collision_rule_thermalized,
        f"CollideSweep{precision_prefix}ThermalizedAVX",
        params_vec
    )

    # generate accessors
    generate_macroscopic_values_accessors(
        ctx,
        collision_rule_unthermalized.method,
        f'macroscopic_values_accessors_{precision_suffix}.h')

    # Boundary conditions
    ubb_dynamic = PatchedUBB(lambda *args: None, dim=3,
                             data_type=data_type_np[ctx.double_accuracy])
    ubb_data_handler = BounceBackSlipVelocityUBB(method.stencil, ubb_dynamic)

    generate_boundary(ctx, f'Dynamic_UBB_{precision_suffix}', ubb_dynamic,
                      method, additional_data_handler=ubb_data_handler,
                      streaming_pattern="push")

    with open(f'Dynamic_UBB_{precision_suffix}.h', 'r+') as f:
        content = f.read()
        f.seek(0)
        f.truncate(0)
        # patch for old versions of pystencils_walberla
        if '#pragma once' not in content:
            content = '#pragma once\n' + content
        # patch for floating point accuracy
        content = content.replace('real_t', data_type_cpp[ctx.double_accuracy])
        f.write(content)

    # communication
    generate_lb_pack_info(
        ctx,
        f'PushPackInfo{precision_prefix}',
        method.stencil,
        fields['pdfs'],
        streaming_pattern='push'
    )

    # Info header containing correct template definitions for stencil and field
    #ctx.write_file(f"InfoHeader{precision_prefix}.h", info_header)

earmark_generated_kernels()
guard_generated_kernels_clang_format()
