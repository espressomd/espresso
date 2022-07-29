#
# Copyright (C) 2020-2022 The ESPResSo project
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
import pystencils_walberla
import pystencils_espresso

import lbmpy
import lbmpy.creationfunctions
import lbmpy.forcemodels
import lbmpy.stencils

import lbmpy_walberla
import lbmpy_espresso

import lees_edwards
import relaxation_rates
import walberla_lbm_generation
import code_generation_context


parser = argparse.ArgumentParser(description='Generate the waLBerla kernels.')
parser.add_argument('--single-precision', action='store_true', required=False,
                    help='Use single-precision')
parser.add_argument("--gpu", action="store_true")
args = parser.parse_args()

if args.gpu:
    target = ps.Target.GPU
else:
    target = ps.Target.CPU

# Make sure we have the correct versions of the required dependencies
SpecifierSet = setuptools.version.pkg_resources.packaging.specifiers.SpecifierSet
for module, requirement in [(ps, '==1.0'), (lbmpy, '==1.0')]:
    assert SpecifierSet(requirement).contains(module.__version__), \
        f"{module.__name__} version {module.__version__} doesn't match requirement {requirement}"


with code_generation_context.CodeGeneration() as ctx:
    ctx.double_accuracy = not args.single_precision
    ctx.cuda = True
    if target == ps.Target.GPU:
        target_prefix = "Cuda"
    else:
        target_prefix = ""

    # vectorization parameters
    cpu_vectorize_info = {
        "instruction_set": "avx",
        "assume_inner_stride_one": True,
        "assume_aligned": True,
        "assume_sufficient_line_padding": False}
    params = {"target": target}
    params_vec = {"target": target, "cpu_vectorize_info": cpu_vectorize_info}

    # codegen configuration
    config = pystencils_espresso.generate_config(ctx, params)

    precision_prefix = pystencils_espresso.precision_prefix[ctx.double_accuracy]
    precision_suffix = pystencils_espresso.precision_suffix[ctx.double_accuracy]
    precision_rng = pystencils_espresso.precision_rng[ctx.double_accuracy]
    kT = sp.symbols('kT')
    stencil = lbmpy.stencils.get_stencil('D3Q19')
    fields = pystencils_espresso.generate_fields(config, stencil)
    force_field = fields['force']

    # LB Method definition
    method = lbmpy.creationfunctions.create_mrt_orthogonal(
        stencil=stencil,
        compressible=True,
        weighted=True,
        relaxation_rates=relaxation_rates.rr_getter,
        force_model=lbmpy.forcemodels.Schiller(force_field.center_vector)
    )

    # generate stream kernels
    pystencils_espresso.generate_stream_sweep(
        ctx,
        method,
        f"{target_prefix}StreamSweep{precision_prefix}",
        params)
    if target == ps.Target.CPU:
        pystencils_espresso.generate_stream_sweep(
            ctx,
            method,
            f"{target_prefix}StreamSweep{precision_prefix}AVX",
            params_vec)

    # generate initial densities
    pdfs_setter = pystencils_espresso.generate_setters(ctx, method, params)
    pystencils_walberla.codegen.generate_sweep(
        ctx,
        f"{target_prefix}InitialPDFsSetter{precision_prefix}",
        pdfs_setter, **params)

    # generate unthermalized collision rule
    collision_rule_unthermalized = lbmpy.creationfunctions.create_lb_collision_rule(
        method,
        zero_centered=True,
        optimization={'cse_global': True,
                      'double_precision': ctx.double_accuracy}
    )
    pystencils_espresso.generate_collision_sweep(
        ctx,
        method,
        collision_rule_unthermalized,
        f"{target_prefix}CollideSweep{precision_prefix}",
        params
    )
    if target == ps.Target.CPU:
        pystencils_espresso.generate_collision_sweep(
            ctx,
            method,
            collision_rule_unthermalized,
            f"{target_prefix}CollideSweep{precision_prefix}AVX",
            params_vec
        )

    # generate unthermalized Lees-Edwards collision rule
    le_config = lbmpy.LBMConfig(stencil=stencil,
                                method=lbmpy.Method.TRT,
                                relaxation_rate=sp.Symbol("omega_shear"),
                                compressible=True,
                                zero_centered=False,
                                force_model=lbmpy.ForceModel.GUO,
                                force=force_field.center_vector,
                                kernel_type='collide_only')
    lbm_opt = lbmpy.LBMOptimisation(symbolic_field=fields["pdfs"])
    le_collision_rule_unthermalized = lbmpy.create_lb_update_rule(
        lbm_config=le_config,
        lbm_optimisation=lbm_opt)

    le_collision_rule_unthermalized = lees_edwards.add_lees_edwards_to_collision(
        config, le_collision_rule_unthermalized,
        fields["pdfs"], stencil, 1)  # shear_dir_normal y
    pystencils_espresso.generate_collision_sweep(
        ctx,
        le_config,
        le_collision_rule_unthermalized,
        f"CollideSweep{precision_prefix}LeesEdwards",
        params
    )
    if target == ps.Target.CPU:
        pystencils_espresso.generate_collision_sweep(
            ctx,
            le_config,
            le_collision_rule_unthermalized,
            f"CollideSweep{precision_prefix}LeesEdwardsAVX",
            params_vec
        )

    # generate thermalized LB
    collision_rule_thermalized = lbmpy.creationfunctions.create_lb_collision_rule(
        method,
        zero_centered=False,
        fluctuating={
            'temperature': kT,
            'block_offsets': 'walberla',
            'rng_node': precision_rng
        },
        optimization={'cse_global': True,
                      'double_precision': ctx.double_accuracy}
    )
    pystencils_espresso.generate_collision_sweep(
        ctx,
        method,
        collision_rule_thermalized,
        f"{target_prefix}CollideSweep{precision_prefix}Thermalized",
        params
    )
    if target == ps.Target.CPU:
        pystencils_espresso.generate_collision_sweep(
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
            config,
            collision_rule_unthermalized.method,
            f'{target_prefix}macroscopic_values_accessors_{precision_suffix}.h')

    # Boundary conditions
    ubb_dynamic = lbmpy_espresso.UBB(
        lambda *args: None, dim=3, data_type=config.data_type.default_factory())
    ubb_data_handler = lbmpy_espresso.BounceBackSlipVelocityUBB(
        method.stencil, ubb_dynamic)

    lbmpy_walberla.generate_boundary(
        ctx, f'{target_prefix}Dynamic_UBB_{precision_suffix}', ubb_dynamic,
        method, additional_data_handler=ubb_data_handler,
        streaming_pattern="push", target=target)

    with open(f'{target_prefix}Dynamic_UBB_{precision_suffix}.h', 'r+') as f:
        content = f.read()
        f.seek(0)
        f.truncate(0)
        # patch for floating point accuracy
        content = content.replace('real_t',
                                  config.data_type.default_factory().c_name)
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
