#
# Copyright (C) 2022-2026 The ESPResSo project
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
import pystencils as ps
import pystencils_walberla
import sympy as sp
import lbmpy
import argparse
import packaging.specifiers
import numpy as np

import pystencils_espresso
import code_generation_context

import ekin
import custom_additional_extensions
import walberla_ek_generation


kernel_codes = "diffusion advection continuity friction_coupling boundary reactions accessors".split()
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

double_precision: bool = not args.single_precision

data_type_cpp = "double" if double_precision else "float"
data_type_np = "float64" if double_precision else "float32"
precision_suffix = pystencils_espresso.precision_suffix[double_precision]
precision_rng = pystencils_espresso.precision_rng_modulo[double_precision]
np2cpp_t = pystencils_espresso.numpy_types_to_cpp_types


def patch_openmp_kernels(content):
    # surrounds omp pragmas with ifdefs
    content = re.sub("^( *#pragma omp .*)$",
                     r"#ifdef _OPENMP\n\1\n#endif", content, flags=re.MULTILINE)
    return content


def patch_unused_direction_arrays_kernel(content, variables):
    for name in variables:
        content = walberla_ek_generation.remove_intermediate_variable(
            content, name)
    return content


def patch_reaction_indexed_kernel(content: str, target_suffix) -> str:
    # replace getData with uncheckedFastGetData
    access_slow = "block->getData<IndexVectors>(indexVectorID);"
    access_fast = "block->uncheckedFastGetData<IndexVectors>(indexVectorID);"
    assert access_slow in content
    content = content.replace(access_slow, access_fast)
    if target_suffix in ["_CUDA"]:
        # replace preprocessor macros and pragmas
        push, pop = custom_additional_extensions.generate_device_preprocessor(
            "reactions", defines=("RESTRICT",))
        content = re.sub(r"#ifdef __GNUC__[\s\S]+?#endif(?=\n\n|\n//)", "", content)  # nopep8
        content = re.sub(r"#ifdef __CUDACC__[\s\S]+?#endif(?=\n\n|\n//)", push, content, 1)  # nopep8
        content = re.sub(r"#ifdef __CUDACC__[\s\S]+?#endif(?=\n\n|\n//)", pop, content, 1)  # nopep8
        assert push in content
        assert pop in content
        # remove dummy assignment
        token = "const int32_t dummy = *((int32_t *  )(& _data_indexVector_10[12*blockDim.x*blockIdx.x + 12*threadIdx.x]));"
        assert token in content
        content = content.replace(token, "")
    else:
        # remove dummy assignment
        token = "const int32_t dummy = *((int32_t *  )(& _data_indexVector[12*ctr_0]));"
        assert token in content
        content = content.replace(token, "")
    content = patch_unused_direction_arrays_kernel(
        content, "cx cy cz invdir".split())
    return content


def patch_dirichlet_boundary_kernel(content: str, target_suffix) -> str:
    content = patch_unused_direction_arrays_kernel(content, ["dir"])
    if target_suffix in ["_CUDA"]:
        # remove unused assignment
        token = "uint8_t * RESTRICT _data_indexVector_112 = _data_indexVector + 12;\n"
        assert token in content
        content = content.replace(token, "")
    return content


def patch_diffusive_flux_elec_kernel(content):
    token = "BlockDataID phiID;\n"
    assert token in content
    content = content.replace(token, "BlockDataID phiID;\npublic:\n  inline void setPhiID(BlockDataID phiID_) { phiID = phiID_; }\nprivate:\n")  # nopep8
    return content


def get_ext_header(target_suffix):
    return {"_CUDA": "h"}.get(target_suffix, "h")


def get_ext_source(target_suffix):
    return {"_CUDA": "cu"}.get(target_suffix, "cpp")


dim: int = 3
target: ps.enums.Target = ps.enums.Target.CPU
if args.gpu:
    target = ps.enums.Target.GPU
if args.kernels == "all" or args.kernels == ["all"]:
    args.kernels = kernel_codes
flux_count: int = 3 ** dim // 2

diffusion = ps.TypedSymbol("D", data_type_np)
kT = ps.TypedSymbol("kT", data_type_np)
valency = ps.TypedSymbol("z", data_type_np)
lb_density = ps.TypedSymbol("rho_lb", data_type_np)
ext_efield = [ps.TypedSymbol(f"f_ext_{i}", data_type_np) for i in range(dim)]

density_field, potential_field, velocity_field, force_field = ps.fields(
    f"rho, phi, u(#), f(#): {data_type_np}[#D]".replace("#", str(dim)), layout='zyxf')
flux_field = ps.fields(
    f"j({flux_count}): {data_type_np}[{dim}D]",
    layout='zyxf',
    field_type=ps.FieldType.STAGGERED_FLUX)

ek = ekin.EK(
    dim=dim,
    density_field=density_field,
    flux_field=flux_field,
    diffusion=diffusion,
    kT=kT,
    lb_density=lb_density,
    velocity_field=velocity_field,
    force_field=force_field,
    potential_field=None,
    valency=None,
    ext_efield=None)
ek_electrostatic = ekin.EK(
    dim=dim,
    density_field=density_field,
    flux_field=flux_field,
    diffusion=diffusion,
    kT=kT,
    velocity_field=velocity_field,
    force_field=force_field,
    potential_field=potential_field,
    valency=valency,
    ext_efield=sp.Matrix(ext_efield))

max_num_reactants: int = 5

react_rhos, orders, stoechom_coefs = [], [], []
for i in range(max_num_reactants):
    react_rhos.append(
        ps.fields(f"rho_{i}: {data_type_np}[#D]".replace("#", str(dim)),
                  layout="zyxf"))
    orders.append(ps.TypedSymbol(f"order_{i}", data_type_np))
    stoechom_coefs.append(ps.TypedSymbol(f"stoech_{i}", data_type_np))
rate_coef = sp.Symbol("rate_coefficient")

reaction_obj = ekin.Reaction(
    species=react_rhos,
    orders=orders,
    stoechom_coefs=stoechom_coefs,
    rate_coef=rate_coef,
)

block_offsets = tuple(
    ps.TypedSymbol(f"block_offset_{i}", np.uint32)
    for i in range(3))

if args.gpu:
    params = {
        "target": target
    }
    cpu_vectorize_info = {}  # dummy handle
    processor_suffix = "_CUDA"
    file_suffix = "cu"
else:
    params = {
        "target": target,
        "cpu_openmp": True,
        "cpu_vectorize_info": {
            "assume_inner_stride_one": False,
        },
    }
    cpu_vectorize_info = params["cpu_vectorize_info"]  # handle to mutable dict
    processor_suffix = ""
    file_suffix = "cpp"


with code_generation_context.CodeGeneration() as ctx:
    ctx.double_accuracy = double_precision
    if target == ps.Target.CPU:
        ctx.openmp = True
    if target == ps.Target.GPU:
        ctx.gpu = True
        ctx.cuda = True

    # codegen configuration
    config = pystencils_espresso.generate_config(ctx, params)

    data_type = "float64" if ctx.double_accuracy else "float32"

    dirichlet_stencil = lbmpy.stencils.LBStencil(
        stencil=((0, 0, 0),), theta0=0.)
    dirichlet = custom_additional_extensions.Dirichlet_Custom(
        lambda *args: None, data_type=data_type_np)
    dirichlet_additional_data = custom_additional_extensions.DirichletAdditionalDataHandler(
        dirichlet_stencil, dirichlet)

    if "diffusion" in args.kernels:
        for midfix, fluctuation in (("", False), ("Thermalized", True)):
            cpu_vectorize_info["cpu_prepend_opt_remove_conditionals"] = False
            class_name = f"DiffusiveFluxKernel{midfix}_{precision_suffix}{processor_suffix}"  # nopep8
            pystencils_walberla.generate_sweep(
                ctx, class_name,
                ek.flux(include_vof=False, include_fluctuations=fluctuation,
                        rng_node=precision_rng),
                staggered=True,
                block_offset=block_offsets if fluctuation else None,
                **params)
            ctx.patch_file(class_name, get_ext_source(
                processor_suffix), patch_openmp_kernels)
            cpu_vectorize_info["cpu_prepend_opt_remove_conditionals"] = False
            class_name = f"DiffusiveFluxKernelWithElectrostatic{midfix}_{precision_suffix}{processor_suffix}"  # nopep8
            pystencils_walberla.generate_sweep(
                ctx, class_name,
                ek_electrostatic.flux(include_vof=False, include_fluctuations=fluctuation,
                                      rng_node=precision_rng),
                staggered=True,
                block_offset=block_offsets if fluctuation else None,
                **params)
            ctx.patch_file(class_name, "h", patch_diffusive_flux_elec_kernel)
            ctx.patch_file(class_name, get_ext_source(
                processor_suffix), patch_openmp_kernels)

    if "advection" in args.kernels:
        def patch_advection_kernel(content, target_suffix):
            if target_suffix in ["_CUDA"]:
                # replace preprocessor macros and pragmas
                token = "#define FUNC_PREFIX __global__"
                assert token in content
                push, _ = custom_additional_extensions.generate_device_preprocessor(
                    "advection", defines=("RESTRICT",))
                content = content.replace(token, f"{token}\n{push}")
            return content

        flux_advection = ps.AssignmentCollection(ek.flux_advection())
        cpu_vectorize_info["cpu_prepend_opt_remove_conditionals"] = False
        class_name = f"AdvectiveFluxKernel_{precision_suffix}{processor_suffix}"  # nopep8
        pystencils_walberla.generate_sweep(
            ctx,
            class_name,
            flux_advection,
            staggered=True,
            **params)
        ctx.patch_file(class_name, get_ext_source(processor_suffix),
                       patch_advection_kernel, processor_suffix)
        ctx.patch_file(class_name, get_ext_source(
            processor_suffix), patch_openmp_kernels)

    if "continuity" in args.kernels:
        class_name = f"ContinuityKernel_{precision_suffix}{processor_suffix}"
        pystencils_walberla.generate_sweep(
            ctx,
            class_name,
            ek.continuity(),
            **params)
        ctx.patch_file(class_name, get_ext_source(
            processor_suffix), patch_openmp_kernels)
    if "friction_coupling" in args.kernels:
        class_name = f"FrictionCouplingKernel_{
            precision_suffix}{processor_suffix}"
        pystencils_walberla.generate_sweep(
            ctx,
            class_name,
            ek.friction_coupling(),
            **params)
        ctx.patch_file(class_name, get_ext_source(
            processor_suffix), patch_openmp_kernels)

    if "boundary" in args.kernels:
        # pylint: disable=unused-argument
        def patch_boundary_header(content, processor_suffix):
            # replace real_t by actual floating-point type
            return content.replace("real_t", f"{np2cpp_t[data_type]}") \
                          .replace("real_c", f"{np2cpp_t[data_type]}_c")

        def patch_boundary_kernel(content, processor_suffix):
            if processor_suffix in ["_CUDA"]:
                # replace preprocessor macros and pragmas
                push, pop = custom_additional_extensions.generate_device_preprocessor(
                    "ubb_boundary", defines=("RESTRICT",))
                content = re.sub(r"#ifdef __GNUC__[\s\S]+?#endif(?=\n\n|\n//)", "", content)  # nopep8
                content = re.sub(r"#ifdef __CUDACC__[\s\S]+?#endif(?=\n\n|\n//)", push, content, 1)  # nopep8
                content = re.sub(r"#ifdef __CUDACC__[\s\S]+?#endif(?=\n\n|\n//)", pop, content, 1)  # nopep8
                assert push in content
                assert pop in content
            content = patch_unused_direction_arrays_kernel(
                content, "cx cy cz invdir".split())
            return content

        # generate dynamic fixed flux
        stencil = lbmpy.LBStencil(stencil="D3Q27")
        dynamic_flux = custom_additional_extensions.Flux(
            stencil, lambda *args: None, dim=3, data_type=data_type_np)
        dynamic_flux_additional_data = custom_additional_extensions.FluxAdditionalDataHandler(
            stencil=stencil, boundary_object=dynamic_flux)
        class_name = f"FixedFlux_{precision_suffix}{processor_suffix}"

        pystencils_walberla.boundary.generate_staggered_flux_boundary(
            generation_context=ctx,
            class_name=class_name,
            boundary_object=dynamic_flux,
            dim=dim,
            neighbor_stencil=stencil,
            index_shape=flux_field.index_shape,
            target=target,
            additional_data_handler=dynamic_flux_additional_data)

        ctx.patch_file(class_name, get_ext_header(processor_suffix),
                       patch_boundary_header, processor_suffix)
        ctx.patch_file(class_name, get_ext_source(processor_suffix),
                       patch_boundary_kernel, processor_suffix)
        ctx.patch_file(class_name, get_ext_source(
            processor_suffix), patch_openmp_kernels)

        # generate dynamic fixed density
        class_name = f"Dirichlet_{precision_suffix}{processor_suffix}"
        pystencils_walberla.boundary.generate_boundary(
            generation_context=ctx,
            class_name=class_name,
            boundary_object=dirichlet,
            additional_data_handler=dirichlet_additional_data,
            field_name="field",
            neighbor_stencil=stencil,
            index_shape=density_field.index_shape,
            target=target)

        ctx.patch_file(class_name, get_ext_header(processor_suffix),
                       patch_boundary_header, processor_suffix)
        ctx.patch_file(class_name, get_ext_source(processor_suffix),
                       patch_boundary_kernel, processor_suffix)
        ctx.patch_file(class_name, get_ext_source(processor_suffix),
                       patch_dirichlet_boundary_kernel, processor_suffix)
        ctx.patch_file(class_name, get_ext_source(
            processor_suffix), patch_openmp_kernels)

    if "reactions" in args.kernels:
        # ek reactions
        for i in range(1, max_num_reactants + 1):
            assignments = list(reaction_obj.generate_reaction(num_reactants=i))
            class_name = f"ReactionKernelBulk_{i}_{precision_suffix}{processor_suffix}"  # nopep8
            pystencils_walberla.generate_sweep(
                generation_context=ctx,
                class_name=class_name,
                target=target,
                assignments=assignments)
            ctx.patch_file(class_name, get_ext_source(
                processor_suffix), patch_openmp_kernels)

            class_name = f"ReactionKernelIndexed_{i}_{precision_suffix}{processor_suffix}"  # nopep8
            custom_additional_extensions.generate_boundary(
                generation_context=ctx,
                stencil=dirichlet_stencil,
                class_name=class_name,
                context_params={
                    "stencil_info": dirichlet_additional_data.stencil_info,
                    "parameters_to_ignore": ["indexVectorSize", "indexVectorID"],
                },
                target=target,
                assignment=assignments,
                template_file="templates/Boundary.tmpl.h")
            ctx.patch_file(class_name, file_suffix,
                           patch_reaction_indexed_kernel, processor_suffix)
            ctx.patch_file(class_name, get_ext_source(
                processor_suffix), patch_openmp_kernels)

        # ek reactions helper functions
        custom_additional_extensions.generate_kernel_selector(
            generation_context=ctx,
            class_name="ReactionKernelBulk",
            max_num_reactants=max_num_reactants,
            precision_suffix=pystencils_espresso.precision_suffix,
            processor_suffix=processor_suffix)
        custom_additional_extensions.generate_kernel_selector(
            generation_context=ctx,
            class_name="ReactionKernelIndexed",
            max_num_reactants=max_num_reactants,
            precision_suffix=pystencils_espresso.precision_suffix,
            processor_suffix=processor_suffix)

    if "accessors" in args.kernels:
        # field accessors
        precision_prefix = pystencils_espresso.precision_prefix[ctx.double_accuracy]
        kernel_name = f"EK_FieldAccessors_{precision_suffix}{processor_suffix}"
        if target == ps.Target.GPU:
            templates = {
                f"{kernel_name}.cuh": "templates/EK_FieldAccessors.tmpl.cuh",
                f"{kernel_name}.cu": "templates/EK_FieldAccessors.tmpl.cu",
            }
        else:
            templates = {
                f"{kernel_name}.h": "templates/EK_FieldAccessors.tmpl.h",
            }
        walberla_ek_generation.generate_accessors(
            ctx, config, templates
        )
