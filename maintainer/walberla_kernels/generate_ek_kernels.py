#
# Copyright (C) 2022-2023 The ESPResSo project
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

import pystencils as ps
import pystencils_walberla
import sympy as sp
import lbmpy
import argparse

import pystencils_espresso
import code_generation_context

import ekin
import custom_additional_extensions


parser = argparse.ArgumentParser(description='Generate the waLBerla kernels.')
parser.add_argument('--single-precision', action='store_true', required=False,
                    help='Use single-precision')
args = parser.parse_args()

double_precision: bool = not args.single_precision

data_type_cpp = "double" if double_precision else "float"
data_type_np = pystencils_espresso.data_type_np[data_type_cpp]
precision_suffix = pystencils_espresso.precision_suffix[double_precision]
precision_rng = pystencils_espresso.precision_rng[double_precision]


def replace_getData_with_uncheckedFastGetData(filename: str) -> None:
    with open(filename, "r+") as f:
        content = f.read()
        f.seek(0)
        f.truncate(0)
        content = content.replace("block->getData<IndexVectors>(indexVectorID);",
                                  "block->uncheckedFastGetData<IndexVectors>(indexVectorID);")
        f.write(content)


dim: int = 3
target: ps.enums.Target = ps.enums.Target.CPU
flux_count: int = 3 ** dim // 2

diffusion = ps.TypedSymbol("D", data_type_np)
kT = ps.TypedSymbol("kT", data_type_np)
valency = ps.TypedSymbol("z", data_type_np)
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

params = {
    "target": target,
    "cpu_vectorize_info": {"assume_inner_stride_one": False}}

with code_generation_context.CodeGeneration() as ctx:
    ctx.double_accuracy = double_precision

    # codegen configuration
    config = pystencils_espresso.generate_config(ctx, params)

    pystencils_walberla.generate_sweep(
        ctx,
        f"DiffusiveFluxKernel_{precision_suffix}",
        ek.flux(include_vof=False, include_fluctuations=False,
                rng_node=precision_rng),
        staggered=True,
        **params)
    pystencils_walberla.generate_sweep(
        ctx,
        f"DiffusiveFluxKernelWithElectrostatic_{precision_suffix}",
        ek_electrostatic.flux(include_vof=False, include_fluctuations=False,
                              rng_node=precision_rng),
        staggered=True,
        **params)
    pystencils_walberla.generate_sweep(
        ctx,
        f"AdvectiveFluxKernel_{precision_suffix}",
        ek.flux_advection(),
        staggered=True,
        **params)
    pystencils_walberla.generate_sweep(
        ctx,
        f"ContinuityKernel_{precision_suffix}",
        ek.continuity(),
        **params)

    pystencils_walberla.generate_sweep(
        ctx,
        f"FrictionCouplingKernel_{precision_suffix}",
        ek.friction_coupling(),
        **params)

    # generate dynamic fixed flux
    stencil = lbmpy.LBStencil(stencil="D3Q27")
    dynamic_flux = custom_additional_extensions.Flux(
        stencil, lambda *args: None, dim=3, data_type=data_type_np)
    dynamic_flux_additional_data = custom_additional_extensions.FluxAdditionalDataHandler(
        stencil=stencil, boundary_object=dynamic_flux)

    pystencils_walberla.generate_staggered_flux_boundary(
        generation_context=ctx,
        class_name=f"FixedFlux_{precision_suffix}",
        boundary_object=dynamic_flux,
        dim=dim,
        neighbor_stencil=stencil,
        index_shape=flux_field.index_shape,
        target=target,
        additional_data_handler=dynamic_flux_additional_data)

    # generate dynamic fixed density
    dirichlet_stencil = lbmpy.stencils.LBStencil(stencil=((0, 0, 0),))
    dirichlet = custom_additional_extensions.Dirichlet_Custom(
        lambda *args: None, data_type=data_type_np)
    dirichlet_additional_data = custom_additional_extensions.DirichletAdditionalDataHandler(
        dirichlet_stencil, dirichlet)

    pystencils_walberla.boundary.generate_boundary(
        generation_context=ctx,
        class_name=f"Dirichlet_{precision_suffix}",
        boundary_object=dirichlet,
        additional_data_handler=dirichlet_additional_data,
        field_name="field",
        neighbor_stencil=stencil,
        index_shape=density_field.index_shape,
        target=target)

    pystencils_walberla.generate_pack_info_from_kernel(
        ctx,
        f"DensityPackInfo_{precision_suffix}",
        ek_electrostatic.continuity(),
        target=target)

    # ek reactions
    for i in range(1, max_num_reactants + 1):
        assignments = list(reaction_obj.generate_reaction(num_reactants=i))
        filename_stem: str = f"ReactionKernelBulk_{i}_{precision_suffix}"
        pystencils_walberla.generate_sweep(
            ctx,
            filename_stem,
            assignments)

        filename_stem: str = f"ReactionKernelIndexed_{i}_{precision_suffix}"
        custom_additional_extensions.generate_boundary(
            generation_context=ctx,
            stencil=dirichlet_stencil,
            class_name=filename_stem,
            dim=dim,
            target=target,
            assignment=assignments)
        replace_getData_with_uncheckedFastGetData(
            filename=f"{filename_stem}.cpp")

    # ek reactions helper functions
    custom_additional_extensions.generate_kernel_selector(
        generation_context=ctx,
        class_name="ReactionKernelBulk",
        max_num_reactants=max_num_reactants,
        precision_suffix=pystencils_espresso.precision_suffix)
    custom_additional_extensions.generate_kernel_selector(
        generation_context=ctx,
        class_name="ReactionKernelIndexed",
        max_num_reactants=max_num_reactants,
        precision_suffix=pystencils_espresso.precision_suffix)
