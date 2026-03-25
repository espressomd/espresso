#
# Copyright (C) 2020-2026 The ESPResSo project
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

import re
import argparse
import packaging.specifiers

import sympy as sp
import numpy as np

import pystencils as ps

import lbmpy

kernel_codes = "packinfo boundary stream_collide init vel_update accessors".split()
parser = argparse.ArgumentParser(description="Generate the waLBerla kernels.")
parser.add_argument("--single-precision", action="store_true", required=False,
                    help="Use single-precision")
parser.add_argument("--gpu", action="store_true")
parser.add_argument("--kernels", nargs="+", type=str, default="all",
                    choices=["all"] + kernel_codes,
                    help="Which kernels to generate")
args = parser.parse_args()

# Make sure we have the correct versions of the required dependencies
for module, requirement in [(ps, "==1.4.0"), (lbmpy, "==1.4.0")]:
    assert packaging.specifiers.SpecifierSet(requirement).contains(module.__version__), \
        f"{module.__name__} version {module.__version__} " \
        f"doesn't match requirement {requirement}"

import pystencils_walberla
import pystencils_espresso
import lbmpy_walberla.additional_data_handler
import lbmpy.creationfunctions
import lbmpy.macroscopic_value_kernels
import lbmpy.forcemodels
import lbmpy.stencils
import lbmpy.enums

import lees_edwards
import relaxation_rates
import walberla_lbm_generation
import code_generation_context
import custom_additional_extensions

if args.gpu:
    target = ps.Target.GPU
else:
    target = ps.Target.CPU
if args.kernels == "all" or args.kernels == ["all"]:
    args.kernels = kernel_codes

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
    parameters["CPU"] = ({"target": target,
                          "cpu_openmp": True}, "")
    parameters["CPU_linear"] = ({"target": target}, "")
    parameters["AVX"] = ({"target": target,
                          "cpu_openmp": True,
                          "cpu_vectorize_info": cpu_vectorize_info}, "AVX")

# global parameters
stencil = lbmpy.stencils.LBStencil(lbmpy.enums.Stencil.D3Q19)
kT = sp.symbols("kT")
lbm_config_kwargs = dict(
    weighted=True,
    compressible=True,
    delta_equilibrium=True,
    zero_centered=True,
)
np2cpp_t = pystencils_espresso.numpy_types_to_cpp_types


class BounceBackSlipVelocityUBB(
        lbmpy_walberla.additional_data_handler.UBBAdditionalDataHandler):
    '''
    Dynamic UBB that implements the bounce-back method with slip velocity.
    '''

    def __init__(self, stencil, boundary_object):
        super().__init__(stencil, boundary_object)
        self.Q = stencil.Q
        self.neighbor_directions = [
            np.array2string(x, separator=",") for x in np.array(
                stencil.stencil_entries).transpose()]


def paramlist(parameters, keys):
    for key in keys:
        if key in parameters:
            yield parameters[key]


def get_ext_header(target_suffix):
    return {"CUDA": "h"}.get(target_suffix, "h")


def get_ext_source(target_suffix):
    return {"CUDA": "cu"}.get(target_suffix, "cpp")


def patch_openmp_kernels(content):
    # surrounds omp pragmas with ifdefs
    content = re.sub("^( *#pragma omp .*)$",
                     r"#ifdef _OPENMP\n\1\n#endif", content, flags=re.MULTILINE)
    return content


def generate_init_kernels(ctx, method):
    precision_prefix = pystencils_espresso.precision_prefix[ctx.double_accuracy]
    for params, target_suffix in paramlist(parameters, (default_key,)):
        stem = f"InitialPDFsSetter{precision_prefix}{target_suffix}"
        pystencils_walberla.generate_sweep(
            ctx,
            stem,
            pystencils_espresso.generate_setters(method, data_type),
            **params)
        ctx.patch_file(stem, get_ext_source(target_suffix),
                       patch_openmp_kernels)


def generate_stream_collide_kernels(ctx, method, data_type, fields):
    precision_prefix = pystencils_espresso.precision_prefix[ctx.double_accuracy]
    optimization = {"cse_global": True,
                    "double_precision": ctx.double_accuracy}
    lbm_config = lbmpy.LBMConfig(stencil=stencil,
                                 streaming_pattern="pull",
                                 **lbm_config_kwargs)
    # Thermalized Stream-Collide Kernels
    precision_rng = pystencils_espresso.precision_rng[ctx.double_accuracy]
    block_offsets = tuple(
        ps.TypedSymbol(f"block_offset_{i}", np.uint32)
        for i in range(3))
    lb_collision_rule_thermalized = lbmpy.creationfunctions.create_lb_collision_rule(
        method,
        lbm_config=lbm_config,
        fluctuating={
            "temperature": kT,
            "block_offsets": block_offsets,
            "rng_node": precision_rng
        },
        optimization=optimization
    )
    for params, target_suffix in paramlist(parameters, ("GPU", "CPU", "AVX")):
        stem = f"StreamCollideSweepThermalized{precision_prefix}{target_suffix}"  # nopep8
        pystencils_espresso.generate_stream_collision_sweep(
            ctx,
            method,
            lbm_config,
            data_type,
            lb_collision_rule_thermalized,
            stem,
            optimization,
            params,
            block_offset=block_offsets,
        )
        ctx.patch_file(stem, get_ext_source(target_suffix),
                       patch_openmp_kernels)
    # Unthermalized Lees-Edwards Stream-Collide Kernels
    shear_dir_normal = 1  # y-axis
    for params, target_suffix in paramlist(parameters, ("GPU", "CPU", "AVX")):
        stem = f"StreamCollideSweepLeesEdwards{precision_prefix}{target_suffix}"  # nopep8
        simd = params.get("cpu_vectorize_info", {}).get("instruction_set")
        le_update_rule_unthermalized = lbmpy.create_lb_collision_rule(
            method,
            lbm_config=lbm_config,
            optimization=optimization
        )
        le_collision_rule_unthermalized = lees_edwards.add_lees_edwards_to_collision(
            config, le_update_rule_unthermalized, fields["pdfs"], stencil,
            shear_dir_normal, simd is not None, True)
        pystencils_espresso.generate_stream_collision_sweep(
            ctx,
            method,
            lbm_config,
            data_type,
            le_collision_rule_unthermalized,
            stem,
            optimization,
            params
        )
        ctx.patch_file(stem, get_ext_source(target_suffix),
                       patch_openmp_kernels)


def generate_vel_update_kernels(ctx, method, fields):
    precision_prefix = pystencils_espresso.precision_prefix[ctx.double_accuracy]
    assignments = lbmpy.macroscopic_value_kernels.macroscopic_values_getter(
        method, None, fields["velocity"], fields["pdfs"],
        use_pre_collision_pdfs=True)
    for params, target_suffix in paramlist(parameters, ("GPU", "CPU", "AVX")):
        stem = f"UpdateVelFromPDF{precision_prefix}{target_suffix}"
        pystencils_walberla.generate_sweep(
            ctx,
            stem,
            assignments,
            **params)
        ctx.patch_file(stem, get_ext_source(target_suffix),
                       patch_openmp_kernels)


def generate_accessors_kernels(ctx, method):
    precision_prefix = pystencils_espresso.precision_prefix[ctx.double_accuracy]
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


def generate_packinfo_kernels(ctx, data_type, fields):
    precision_prefix = pystencils_espresso.precision_prefix[ctx.double_accuracy]
    assignments = pystencils_espresso.generate_pack_info_pdfs_field_assignments(
        fields, streaming_pattern="pull")
    spec = pystencils_espresso.generate_pack_info_field_specifications(
        stencil, data_type, fields["force"].layout)

    def patch_packinfo_header(content, target_suffix):
        if target_suffix in ["", "AVX"]:
            # remove todo comment
            token = "\n       //TODO: optimize by generating kernel for this case\n"
            assert token in content
            content = content.replace(token, "\n")
        elif target_suffix in ["CUDA"]:
            # replace preprocessor macros and pragmas
            token = "#define FUNC_PREFIX __global__"
            assert token in content
            content = content.replace(token, "")
            content = re.sub(r"#ifdef __GNUC__[\s\S]+?#endif\n\n", "", content)
        return content

    def patch_packinfo_kernel(content, target_suffix):
        if target_suffix in ["CUDA"]:
            # replace preprocessor macros and pragmas
            token = "#define FUNC_PREFIX __global__"
            assert token in content
            push, _ = custom_additional_extensions.generate_device_preprocessor(
                "packinfo", defines=("RESTRICT",))
            content = content.replace(token, f"{token}\n{push}")
            # add missing includes
            token = '#include "PackInfo'
            assert token in content
            content = content.replace(token, f'#include "core/DataTypes.h"\n#include "core/cell/CellInterval.h"\n#include "domain_decomposition/IBlock.h"\n#include "stencil/Directions.h"\n\n{token}')  # nopep8
        return content

    for params, target_suffix in paramlist(parameters, ["CPU_linear", "GPU"]):
        pystencils_walberla.generate_pack_info_from_kernel(
            ctx, f"PackInfoPdf{precision_prefix}{target_suffix}", assignments,
            kind="pull", **params)
        pystencils_walberla.generate_pack_info(
            ctx, f"PackInfoVec{precision_prefix}{target_suffix}", spec, **params)
        for suffix in ["Pdf", "Vec"]:
            class_name = f"PackInfo{suffix}{precision_prefix}{target_suffix}"
            ctx.patch_file(class_name, get_ext_header(target_suffix),
                           patch_packinfo_header, target_suffix)
            ctx.patch_file(class_name, get_ext_source(target_suffix),
                           patch_packinfo_kernel, target_suffix)


def generate_boundary_kernels(ctx, method, data_type):
    precision_prefix = pystencils_espresso.precision_prefix[ctx.double_accuracy]
    ubb_dynamic = lbmpy.boundaries.UBB(
        lambda *args: None, dim=3, data_type=data_type, calculate_force_on_boundary=True)
    ubb_data_handler = BounceBackSlipVelocityUBB(method.stencil, ubb_dynamic)

    # pylint: disable=unused-argument
    def patch_boundary_header(content, target_suffix):
        # replace real_t by actual floating-point type
        return content.replace("real_t", f"{np2cpp_t[data_type]}") \
                      .replace("real_c", f"{np2cpp_t[data_type]}_c")

    def patch_boundary_kernel(content, target_suffix):
        if target_suffix in ["CUDA"]:
            # replace preprocessor macros and pragmas
            push, pop = custom_additional_extensions.generate_device_preprocessor(
                "ubb_boundary", defines=("RESTRICT",))
            content = re.sub(r"#ifdef __GNUC__[\s\S]+?#endif(?=\n\n|\n//)", "", content)  # nopep8
            content = re.sub(r"#ifdef __CUDACC__[\s\S]+?#endif(?=\n\n|\n//)", push, content, 1)  # nopep8
            content = re.sub(r"#ifdef __CUDACC__[\s\S]+?#endif(?=\n\n|\n//)", pop, content, 1)  # nopep8
            assert push in content
            assert pop in content
        return content

    for _, target_suffix in paramlist(parameters, ("CPU", "GPU")):
        class_name = f"DynamicUBB{precision_prefix}{target_suffix}"
        custom_additional_extensions.generate_lb_boundary(
            ctx, class_name, ubb_dynamic, method,
            additional_data_handler=ubb_data_handler,
            streaming_pattern="pull", target=target,
            template_file="templates/Boundary.tmpl.h")
        ctx.patch_file(class_name, get_ext_header(target_suffix),
                       patch_boundary_header, target_suffix)
        ctx.patch_file(class_name, get_ext_source(target_suffix),
                       patch_boundary_kernel, target_suffix)
        ctx.patch_file(class_name, get_ext_source(target_suffix),
                       patch_openmp_kernels)


with code_generation_context.CodeGeneration() as ctx:
    ctx.double_accuracy = not args.single_precision
    if target == ps.Target.CPU:
        ctx.openmp = True
    if target == ps.Target.GPU:
        ctx.gpu = True
        ctx.cuda = True

    # codegen configuration
    config = pystencils_espresso.generate_config(
        ctx, parameters[default_key][0])
    data_type = "float64" if ctx.double_accuracy else "float32"
    fields = pystencils_espresso.generate_fields(stencil, data_type)

    # LB Method definition
    method = lbmpy.creationfunctions.create_mrt_orthogonal(
        stencil=stencil,
        relaxation_rates=relaxation_rates.rr_getter,
        force_model=lbmpy.forcemodels.Schiller(fields["force"].center_vector),
        **lbm_config_kwargs
    )

    if "init" in args.kernels:
        generate_init_kernels(ctx, method)
    if "stream_collide" in args.kernels:
        generate_stream_collide_kernels(ctx, method, data_type, fields)
    if "vel_update" in args.kernels:
        generate_vel_update_kernels(ctx, method, fields)
    if "accessors" in args.kernels:
        generate_accessors_kernels(ctx, method)
    if "packinfo" in args.kernels:
        generate_packinfo_kernels(ctx, data_type, fields)
    if "boundary" in args.kernels:
        generate_boundary_kernels(ctx, method, data_type)
