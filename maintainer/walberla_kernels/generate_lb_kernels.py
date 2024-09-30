#
# Copyright (C) 2020-2023 The ESPResSo project
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
import packaging.specifiers

import sympy as sp
import numpy as np

import pystencils as ps
import pystencils_walberla
import pystencils_espresso

import lbmpy
import lbmpy.creationfunctions
import lbmpy.forcemodels
import lbmpy.stencils
import lbmpy.enums

import lbmpy_walberla
import lbmpy_espresso

import lees_edwards
import relaxation_rates
import walberla_lbm_generation
import code_generation_context

parser = argparse.ArgumentParser(description="Generate the waLBerla kernels.")
parser.add_argument("--single-precision", action="store_true", required=False,
                    help="Use single-precision")
parser.add_argument("--gpu", action="store_true")
args = parser.parse_args()

if args.gpu:
    target = ps.Target.GPU
else:
    target = ps.Target.CPU

# Make sure we have the correct versions of the required dependencies
for module, requirement in [(ps, "==1.3.3"), (lbmpy, "==1.3.3")]:
    assert packaging.specifiers.SpecifierSet(requirement).contains(module.__version__), \
        f"{module.__name__} version {module.__version__} " \
        f"doesn't match requirement {requirement}"


def paramlist(parameters, keys):
    for key in keys:
        if key in parameters:
            yield parameters[key]


with code_generation_context.CodeGeneration() as ctx:
    ctx.double_accuracy = not args.single_precision
    if target == ps.Target.GPU:
        ctx.gpu = True
        ctx.cuda = True

    # vectorization parameters
    parameters = {}
    if target == ps.Target.GPU:
        default_key = "GPU"
        parameters["GPU"] = ({"target": target}, "CUDA")
    else:
        default_key = "CPU"
        cpu_vectorize_info = {
            "instruction_set": "avx",
            "assume_inner_stride_one": True,
            "assume_aligned": True,
            "assume_sufficient_line_padding": False}
        parameters["CPU"] = ({"target": target}, "")
        parameters["AVX"] = ({"target": target,
                             "cpu_vectorize_info": cpu_vectorize_info}, "AVX")

    # codegen configuration
    config = pystencils_espresso.generate_config(
        ctx, parameters[default_key][0])

    precision_prefix = pystencils_espresso.precision_prefix[ctx.double_accuracy]
    precision_suffix = pystencils_espresso.precision_suffix[ctx.double_accuracy]
    precision_rng = pystencils_espresso.precision_rng[ctx.double_accuracy]
    kT = sp.symbols("kT")
    stencil = lbmpy.stencils.LBStencil(lbmpy.enums.Stencil.D3Q19)
    fields = pystencils_espresso.generate_fields(config, stencil)
    force_field = fields["force"]
    lbm_opt = lbmpy.LBMOptimisation(symbolic_field=fields["pdfs"])
    streaming_pattern = "push"

    # LB Method definition
    method = lbmpy.creationfunctions.create_mrt_orthogonal(
        stencil=stencil,
        compressible=True,
        weighted=True,
        relaxation_rates=relaxation_rates.rr_getter,
        force_model=lbmpy.forcemodels.Schiller(force_field.center_vector)
    )

    # generate stream kernels
    for params, target_suffix in paramlist(parameters, ("GPU", "CPU", "AVX")):
        pystencils_espresso.generate_stream_sweep(
            ctx,
            method,
            f"StreamSweep{precision_prefix}{target_suffix}",
            params)

    # generate initial densities
    for params, target_suffix in paramlist(parameters, (default_key,)):
        pystencils_walberla.generate_sweep(
            ctx,
            f"InitialPDFsSetter{precision_prefix}{target_suffix}",
            pystencils_espresso.generate_setters(ctx, method, params),
            **params)

    # generate unthermalized Lees-Edwards collision rule
    le_config = lbmpy.LBMConfig(stencil=stencil,
                                method=lbmpy.Method.TRT,
                                relaxation_rate=sp.Symbol("omega_shear"),
                                compressible=True,
                                zero_centered=False,
                                force_model=lbmpy.ForceModel.GUO,
                                force=force_field.center_vector,
                                kernel_type="collide_only")
    le_update_rule_unthermalized = lbmpy.create_lb_update_rule(
        lbm_config=le_config,
        lbm_optimisation=lbm_opt)
    le_collision_rule_unthermalized = lees_edwards.add_lees_edwards_to_collision(
        config, le_update_rule_unthermalized,
        fields["pdfs"], stencil, 1)  # shear_dir_normal y
    for params, target_suffix in paramlist(parameters, ("GPU", "CPU", "AVX")):
        pystencils_espresso.generate_collision_sweep(
            ctx,
            le_config,
            le_collision_rule_unthermalized,
            f"CollideSweep{precision_prefix}LeesEdwards{target_suffix}",
            params
        )

    block_offsets = tuple(
        ps.TypedSymbol(f"block_offset_{i}", np.uint32)
        for i in range(3))

    # generate thermalized LB collision rule
    lb_collision_rule_thermalized = lbmpy.creationfunctions.create_lb_collision_rule(
        method,
        zero_centered=False,
        fluctuating={
            "temperature": kT,
            "block_offsets": block_offsets,
            "rng_node": precision_rng
        },
        optimization={"cse_global": True,
                      "double_precision": ctx.double_accuracy}
    )
    for params, target_suffix in paramlist(parameters, ("GPU", "CPU", "AVX")):
        stem = f"CollideSweep{precision_prefix}Thermalized{target_suffix}"
        pystencils_espresso.generate_collision_sweep(
            ctx,
            method,
            lb_collision_rule_thermalized,
            stem,
            params,
            block_offset=block_offsets,
        )

    # generate accessors
    for _, target_suffix in paramlist(parameters, ("GPU", "CPU")):
        stem = f"FieldAccessors{precision_prefix}{target_suffix}"
        if target == ps.Target.GPU:
            templates = {
                f"{stem}.cuh": "templates/FieldAccessors.tmpl.cuh",
                f"{stem}.cu": "templates/FieldAccessors.tmpl.cu",
            }
        else:
            templates = {
                f"{stem}.h": "templates/FieldAccessors.tmpl.h",
            }
        walberla_lbm_generation.generate_macroscopic_values_accessors(
            ctx, config, method, templates
        )

    # generate PackInfo
    assignments = pystencils_espresso.generate_pack_info_pdfs_field_assignments(
        fields, streaming_pattern="pull")
    spec = pystencils_espresso.generate_pack_info_vector_field_specifications(
        config, stencil, force_field.layout)
    for params, target_suffix in paramlist(parameters, ["CPU"]):
        pystencils_walberla.generate_pack_info_from_kernel(
            ctx, f"PackInfoPdf{precision_prefix}{target_suffix}", assignments,
            kind="pull", **params)
        pystencils_walberla.generate_pack_info(
            ctx, f"PackInfoVec{precision_prefix}{target_suffix}", spec, **params)
        if target_suffix == "CUDA":
            continue
        token = "\n       //TODO: optimize by generating kernel for this case\n"
        for field_suffix in ["Pdf", "Vec"]:
            class_name = f"PackInfo{field_suffix}{precision_prefix}{target_suffix}"  # nopep8
            with open(f"{class_name}.h", "r+") as f:
                content = f.read()
                assert token in content
                content = content.replace(token, "\n")
                f.seek(0)
                f.truncate()
                f.write(content)

    # boundary conditions
    ubb_dynamic = lbmpy_espresso.UBB(
        lambda *args: None, dim=3, data_type=config.data_type.default_factory())
    ubb_data_handler = lbmpy_espresso.BounceBackSlipVelocityUBB(
        method.stencil, ubb_dynamic)

    for _, target_suffix in paramlist(parameters, ("GPU", "CPU")):
        lbmpy_walberla.generate_boundary(
            ctx, f"Dynamic_UBB_{precision_suffix}{target_suffix}", ubb_dynamic,
            method, additional_data_handler=ubb_data_handler,
            streaming_pattern=streaming_pattern, target=target)

        with open(f"Dynamic_UBB_{precision_suffix}{target_suffix}.h", "r+") as f:
            content = f.read()
            f.seek(0)
            f.truncate(0)
            # patch for floating point accuracy
            content = content.replace("real_t",
                                      config.data_type.default_factory().c_name)
            f.write(content)
