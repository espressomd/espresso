#
#  Copyright (C) 2022 The ESPResSo project
#
#  This file is part of ESPResSo.
#
#  ESPResSo is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ESPResSo is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
import numpy as np
import pystencils as ps
import pystencils_walberla
import sympy as sp

import pystencils_espresso

import ekin
import lbmpy

import custom_additional_extensions

import code_generation_context

import argparse

parser = argparse.ArgumentParser(description='Generate the waLBerla kernels.')
parser.add_argument('--single-precision', action='store_true', required=False,
                    help='Use single-precision')
args = parser.parse_args()

double_precision: bool = not args.single_precision

data_type_cpp = pystencils_espresso.data_type_cpp[double_precision]
data_type_np = pystencils_espresso.data_type_np[double_precision]
precision_suffix = pystencils_espresso.precision_suffix[double_precision]
precision_rng = pystencils_espresso.precision_rng[double_precision]


def replace_real_t_and_add_pragma(filename: str, data_type: str) -> None:
    with open(filename, "r+") as f:
        content = f.read()
        f.seek(0)
        f.truncate(0)
        # patch for old versions of pystencils_walberla
        if '#pragma once' not in content:
            content = '#pragma once\n' + content
        # patch for floating point accuracy
        content = content.replace('real_t', data_type)
        content = content.replace('real_c', f"numeric_cast<{data_type}>")
        f.write(content)


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
ext_efield = []
for i in range(dim):
    ext_efield.append(ps.TypedSymbol(f"f_ext_{i}", data_type_np))

rho, phi, u, f = ps.fields(
    f"rho, phi, u(#), f(#): {data_type_np}[#D]".replace(
        "#", str(dim)), layout='zyxf')
j = ps.fields(
    f"j({flux_count}): {data_type_np}[{dim}D]",
    layout='zyxf',
    field_type=ps.FieldType.STAGGERED_FLUX)

ek = ekin.EK(dim, rho, j, diffusion, kT=kT, u=u, f=f)
ek_electrostatic = ekin.EK(dim, rho, j, diffusion, kT=kT, u=u, f=f, phi=phi, valency=valency,
                           ext_efield=sp.Matrix(ext_efield))

max_num_reactants: int = 5

react_rhos, orders, stoechom_coefs = [], [], []
for i in range(max_num_reactants):
    react_rhos.append(
        ps.fields(
            f"rho_{i}: {data_type_np}[#D]".replace(
                "#",
                str(dim)),
            layout="zyxf"))
    orders.append(ps.TypedSymbol(f"order_{i}", data_type_np))
    stoechom_coefs.append(ps.TypedSymbol(f"stoech_{i}", data_type_np))
rate_coef = sp.Symbol("rate_coefficient")

coordinate_names = ("x", "y", "z")[:dim]
coordinate_size = np.dtype(np.int32).itemsize

index_struct_dtype = np.dtype(
    [(name, np.int32) for name in coordinate_names], align=True
)

index_field = ps.Field(
    "indexVector",
    ps.FieldType.INDEXED,
    index_struct_dtype,
    layout=[0],
    shape=(
        ps.data_types.TypedSymbol(
            "indexVectorSize", ps.data_types.create_type(np.int64)
        ),
        1,
    ),
    strides=(1, 1),
)

config = ps.CreateKernelConfig(
    index_fields=[index_field], coordinate_names=coordinate_names
)

reaction_obj = ekin.Reaction(
    species=react_rhos,
    orders=orders,
    stoechom_coefs=stoechom_coefs,
    rate_coef=rate_coef,
)

params = {'cpu_vectorize_info': {'assume_inner_stride_one': False}}

with code_generation_context.CodeGeneration() as ctx:
    ctx.double_accuracy = double_precision

    # this only necessary because of a bug in 0.4.4 that was fixed shortly
    # after..
    bug_params = {"omp_single_loop": False}

    pystencils_walberla.generate_sweep(ctx, f"DiffusiveFluxKernel_{precision_suffix}",
                                       ek.flux(include_vof=False, include_fluctuations=False), staggered=True, **params,
                                       **bug_params)
    pystencils_walberla.generate_sweep(ctx, f"DiffusiveFluxKernelWithElectrostatic_{precision_suffix}",
                                       ek_electrostatic.flux(
                                           include_vof=False, include_fluctuations=False),
                                       staggered=True, **params, **bug_params)
    pystencils_walberla.generate_sweep(ctx, f"AdvectiveFluxKernel_{precision_suffix}", ek.flux_advection(),
                                       staggered=True, **params, **bug_params)
    pystencils_walberla.generate_sweep(
        ctx,
        f"ContinuityKernel_{precision_suffix}",
        ek.continuity(),
        **params)
    pystencils_walberla.generate_sweep(ctx, f"FrictionCouplingKernel_{precision_suffix}", ek.friction_coupling(),
                                       **params)

    # generate dynamic fixed flux
    stencil = lbmpy.LBStencil(stencil="D3Q27")
    dynamic_flux = custom_additional_extensions.Flux(
        stencil, lambda *args: None, dim=3, data_type=data_type_np)
    dynamic_flux_additional_data = custom_additional_extensions.FluxAdditionalDataHandler(
        stencil=stencil, boundary_object=dynamic_flux)

    pystencils_walberla.generate_staggered_flux_boundary(generation_context=ctx,
                                                         class_name=f"FixedFlux_{precision_suffix}",
                                                         boundary_object=dynamic_flux,
                                                         dim=dim,
                                                         neighbor_stencil=stencil,
                                                         index_shape=j.index_shape,
                                                         target=target,
                                                         additional_data_handler=dynamic_flux_additional_data)

    filename_stem: str = f"FixedFlux_{precision_suffix}"
    replace_real_t_and_add_pragma(
        filename=f"{filename_stem}.h",
        data_type=data_type_cpp)

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
        field_name='field',
        neighbor_stencil=stencil,
        index_shape=rho.index_shape,
        target=target)

    filename_stem: str = f"Dirichlet_{precision_suffix}"
    replace_real_t_and_add_pragma(
        filename=f"{filename_stem}.h",
        data_type=data_type_cpp)

    pystencils_walberla.generate_pack_info_from_kernel(ctx, f"DensityPackInfo_{precision_suffix}",
                                                       ek_electrostatic.continuity(),
                                                       target=target)

    # ek reactions
    for i in range(1, max_num_reactants + 1):
        assignments = list(reaction_obj.generate_reaction(num_reactants=i))
        # this one is correct, as expected
        pystencils_walberla.generate_sweep(
            ctx, f"ReactionKernelBulk_{i}_{precision_suffix}", assignments)

        filename_stem: str = f"ReactionKernelIndexed_{i}_{precision_suffix}"
        custom_additional_extensions.generate_boundary(generation_context=ctx, stencil=dirichlet_stencil,
                                                       class_name=filename_stem,
                                                       dim=dim, target=config.target, assignment=assignments)
        replace_getData_with_uncheckedFastGetData(
            filename=f"{filename_stem}.cpp")
