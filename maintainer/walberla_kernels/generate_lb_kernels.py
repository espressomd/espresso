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

import argparse
import setuptools

import sympy as sp
import pystencils as ps
from pystencils_walberla import CodeGeneration, codegen

import lbmpy
import lbmpy.advanced_streaming.indexing
import lbmpy.boundaries
import lbmpy.creationfunctions
import lbmpy.fieldaccess
import lbmpy.forcemodels
import lbmpy.macroscopic_value_kernels
import lbmpy.stencils
import lbmpy.updatekernels

import lbmpy_walberla
import lbmpy_walberla.additional_data_handler

import lees_edwards
import relaxation_rates
import walberla_lbm_generation
import post_process_kernels


parser = argparse.ArgumentParser(description='Generate the waLBerla kernels.')
parser.add_argument('--single-precision', action='store_true', required=False,
                    help='Use single-precision')
parser.add_argument("--gpu", action="store_true")
args = parser.parse_args()

if args.gpu:
    target = ps.Target.GPU
else:
    target = ps.Target.CPU
data_type_cpp = {True: 'double', False: 'float'}
data_type_np = {True: 'float64', False: 'float32'}

# Make sure we have the correct versions of the required dependencies
SpecifierSet = setuptools.version.pkg_resources.packaging.specifiers.SpecifierSet
for module, requirement in [(ps, '==0.4.4'), (lbmpy, '==0.4.4')]:
    assert SpecifierSet(requirement).contains(module.__version__), \
        f"{module.__name__} version {module.__version__} doesn't match requirement {requirement}"


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
    collide_update_rule = lbmpy.updatekernels.create_lbm_kernel(
        collision_rule,
        fields['pdfs'],
        fields['pdfs_tmp'],
        lbmpy.fieldaccess.CollideOnlyInplaceAccessor())
    collide_ast = ps.create_kernel(collide_update_rule, **params,
                                   data_type=data_type_np[ctx.double_accuracy])
    collide_ast.function_name = 'kernel_collide'
    collide_ast.assumed_inner_stride_one = True
    codegen.generate_sweep(ctx, class_name, collide_ast, target=target)


def generate_stream_sweep(ctx, lb_method, class_name, params):
    # Symbols for PDF (twice, due to double buffering)
    fields = generate_fields(ctx, lb_method.stencil)

    # Generate stream kernel
    stream_update_rule = lbmpy.updatekernels.create_stream_pull_with_output_kernel(
        lb_method, fields['pdfs'], fields['pdfs_tmp'],
        output={'velocity': fields['velocity']})
    stream_ast = ps.create_kernel(stream_update_rule, **params,
                                  data_type=data_type_np[ctx.double_accuracy])
    stream_ast.function_name = 'kernel_stream'
    stream_ast.assumed_inner_stride_one = True
    codegen.generate_sweep(
        ctx, class_name, stream_ast, field_swaps=[(fields['pdfs'], fields['pdfs_tmp'])], target=target)


def generate_setters(lb_method):
    fields = generate_fields(ctx, lb_method.stencil)

    initial_rho = sp.Symbol('rho_0')
    pdfs_setter = lbmpy.macroscopic_value_kernels.macroscopic_values_setter(
        lb_method,
        initial_rho,
        fields['velocity'].center_vector,
        fields['pdfs'].center_vector)
    return pdfs_setter



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

        neighbor_offset = lbmpy.advanced_streaming.indexing.NeighbourOffsetArrays.neighbour_offset(
            dir_symbol, lb_method.stencil)

        assignment = assignments[-1]
        assert assignment.lhs.field == f_in
        out.append(ps.Assignment(assignment.lhs.get_shifted(*neighbor_offset),
                                 assignment.rhs - f_out(dir_symbol) + f_in(dir_symbol)))
        return out


post_process_kernels.adapt_pystencils()


class PatchedCodeGeneration(CodeGeneration):

    def __init__(self):
        import sys
        old_sys_argv = sys.argv
        sys.argv = sys.argv[:1]
        super().__init__()
        sys.argv = old_sys_argv


with PatchedCodeGeneration() as ctx:
    ctx.double_accuracy = not args.single_precision
    ctx.cuda = True
    if target == ps.Target.GPU:
        target_prefix = "Cuda"
    else:
        target_prefix = ""
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
    stencil = lbmpy.stencils.get_stencil('D3Q19')
    fields = generate_fields(ctx, stencil)
    force_field = fields['force']

    # Vectorization parameters
    cpu_vectorize_info = {
        "instruction_set": "avx",
        "assume_inner_stride_one": True,
        "assume_aligned": True,
        "assume_sufficient_line_padding": False}
    params = {"target": target}
    params_vec = {"target": target, "cpu_vectorize_info": cpu_vectorize_info}

    # LB Method definition
    method = lbmpy.creationfunctions.create_mrt_orthogonal(
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
        f"{target_prefix}StreamSweep{precision_prefix}",
        params)
    if target == ps.Target.CPU:
        generate_stream_sweep(
            ctx,
            method,
            f"{target_prefix}StreamSweep{precision_prefix}AVX",
            params_vec)

    # generate initial densities
    pdfs_setter = generate_setters(method)
    codegen.generate_sweep(
        ctx,
        f"{target_prefix}InitialPDFsSetter{precision_prefix}",
        pdfs_setter, target=target)

    # generate unthermalized collision rule
    collision_rule_unthermalized = lbmpy.creationfunctions.create_lb_collision_rule(
        method,
        optimization={'cse_global': True,
                      'double_precision': ctx.double_accuracy}
    )
    generate_collision_sweep(
        ctx,
        method,
        collision_rule_unthermalized,
        f"{target_prefix}CollideSweep{precision_prefix}",
        params
    )
    if target == ps.Target.CPU:
        generate_collision_sweep(
            ctx,
            method,
            collision_rule_unthermalized,
            f"{target_prefix}CollideSweep{precision_prefix}AVX",
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
    collision_rule_thermalized = lbmpy.creationfunctions.create_lb_collision_rule(
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
        f"{target_prefix}CollideSweep{precision_prefix}Thermalized",
        params
    )
    if target == ps.Target.CPU:
        generate_collision_sweep(
            ctx,
            method,
            collision_rule_thermalized,
            f"{target_prefix}CollideSweep{precision_prefix}ThermalizedAVX",
            params_vec
        )

    # generate accessors
    if target == ps.Target.CPU:
        walberla_lbm_generation.generate_macroscopic_values_accessors(
            ctx,
            collision_rule_unthermalized.method,
            f'{target_prefix}macroscopic_values_accessors_{precision_suffix}.h')

    # Boundary conditions
    ubb_dynamic = PatchedUBB(lambda *args: None, dim=3,
                             data_type=data_type_np[ctx.double_accuracy])
    ubb_data_handler = BounceBackSlipVelocityUBB(method.stencil, ubb_dynamic)

    lbmpy_walberla.generate_boundary(
        ctx, f'{target_prefix}Dynamic_UBB_{precision_suffix}', ubb_dynamic,
        method, additional_data_handler=ubb_data_handler,
        streaming_pattern="push", target=target)

    with open(f'{target_prefix}Dynamic_UBB_{precision_suffix}.h', 'r+') as f:
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
    lbmpy_walberla.generate_lb_pack_info(
        ctx,
        f'{target_prefix}PushPackInfo{precision_prefix}',
        method.stencil,
        fields['pdfs'],
        streaming_pattern='push',
        target=target
    )

post_process_kernels.earmark_generated_kernels()
post_process_kernels.guard_generated_kernels_clang_format()
