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
import json
import jinja2
import argparse

import sympy as sp
import pystencils as ps
from lbmpy.creationfunctions import create_lb_collision_rule, create_mrt_orthogonal, force_model_from_string
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
parser.add_argument('codegen_cfg', type=str, nargs='?',
                    help='codegen configuration file (JSON format)')
args = parser.parse_args()


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


def generate_fields(ctx, stencil):
    dtype = 'float64' if ctx.double_accuracy else 'float32'
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
                                   data_type='double' if ctx.double_accuracy else 'float')
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
                                  data_type='double' if ctx.double_accuracy else 'float')
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


adapt_pystencils()


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


class PatchedCodeGeneration(CodeGeneration):
    '''
    Code generator class that reads the configuration from a file.
    '''

    def __init__(self, codegen_cfg):
        if codegen_cfg:
            assert sys.argv[1] == codegen_cfg
            del sys.argv[1]
        super().__init__()
        if codegen_cfg:
            cmake_map = {"ON": True, "1": True, "YES": True,
                         "OFF": False, "0": False, "NO": False}
            with open(codegen_cfg) as f:
                cmake_vars = json.loads(f.read())['CMAKE_VARS']
                for key, value in cmake_vars.items():
                    cmake_vars[key] = cmake_map.get(value, value)
            self.context = type(self.context)(cmake_vars)


with PatchedCodeGeneration(args.codegen_cfg) as ctx:
    kT = sp.symbols("kT")
    stencil = get_stencil('D3Q19')
    fields = generate_fields(ctx, stencil)
    force_field = fields['force']

    # Vectorization parameters
    cpu_vectorize_info = {
        "instruction_set": "avx",
        "assume_inner_stride_one": True,
        "assume_aligned": True,
        "assume_sufficient_line_padding": True}
    params = {"target": "cpu"}
    params_vec = {"target": "cpu", "cpu_vectorize_info": cpu_vectorize_info}

    # LB Method definition
    method = create_mrt_orthogonal(
        stencil=stencil,
        compressible=True,
        weighted=True,
        relaxation_rate_getter=relaxation_rates.rr_getter,
        force_model=force_model_from_string(
            'schiller', force_field.center_vector)
    )
    generate_stream_sweep(ctx, method, "StreamSweep", params)
    generate_stream_sweep(ctx, method, "StreamSweepAVX", params_vec)

    # generate initial densities
    pdfs_setter = generate_setters(method)
    codegen.generate_sweep(ctx, "InitialPDFsSetter", pdfs_setter)

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
        "CollideSweep",
        {}
    )
    generate_collision_sweep(
        ctx,
        method,
        collision_rule_unthermalized,
        "CollideSweepAVX",
        {"cpu_vectorize_info": cpu_vectorize_info}
    )

    # generate unthermalized Lees-Edwards collision rule
    u_p = lees_edwards.velocity_shift()
    velocity_field = ps.fields("velocity(3): [3D]", layout='fzyx')
    le_collision_rule_unthermalized = create_lb_collision_rule(
        method,
        velocity_input=velocity_field.center_vector + sp.Matrix([u_p, 0, 0]),
        optimization={'cse_global': True}
    )
    lees_edwards.modify_method(le_collision_rule_unthermalized, 0, 1)
    generate_collision_sweep(
        ctx,
        method,
        le_collision_rule_unthermalized,
        "CollideSweepLeesEdwards",
        {}
    )
    generate_collision_sweep(
        ctx,
        method,
        le_collision_rule_unthermalized,
        "CollideSweepLeesEdwardsAVX",
        {"cpu_vectorize_info": cpu_vectorize_info}
    )

    # generate thermalized LB
    collision_rule_thermalized = create_lb_collision_rule(
        method,
        fluctuating={
            'temperature': kT,
            'block_offsets': 'walberla',
            'rng_node': ps.rng.PhiloxTwoDoubles
        },
        optimization={'cse_global': True,
                      'double_precision': ctx.double_accuracy}
    )
    generate_collision_sweep(
        ctx,
        method,
        collision_rule_thermalized,
        "CollideSweepThermalized",
        params
    )
    generate_collision_sweep(
        ctx,
        method,
        collision_rule_thermalized,
        "CollideSweepThermalizedAVX",
        params_vec
    )

    # generate accessors
    generate_macroscopic_values_accessors(
        ctx,
        collision_rule_unthermalized,
        field_layout="fzyx")

    # Boundary conditions
    ubb_dynamic = PatchedUBB(lambda *args: None, dim=3)
    ubb_data_handler = BounceBackSlipVelocityUBB(method.stencil, ubb_dynamic)

    generate_boundary(ctx, 'Dynamic_UBB', ubb_dynamic, method,
                      additional_data_handler=ubb_data_handler,
                      streaming_pattern="push")

    # patch for old versions of pystencils_walberla
    with open('Dynamic_UBB.h', 'r+') as f:
        content = f.read()
        if '#pragma once' not in content:
            f.seek(0)
            f.write('#pragma once\n' + content)

    # communication
    generate_lb_pack_info(
        ctx,
        'PushPackInfo',
        method.stencil,
        fields['pdfs'],
        streaming_pattern='push'
    )

    # Info header containing correct template definitions for stencil and field
    #ctx.write_file("InfoHeader.h", info_header)

earmark_generated_kernels()
