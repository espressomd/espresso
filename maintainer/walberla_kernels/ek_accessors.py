#
# Copyright (C) 2021-2026 The ESPResSo project
# Copyright (C) 2020-2022 The waLBerla project
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

"""
Catalogue of the field accessors of the electrokinetics of the waLBerla
bridge (``walberla::ek::accessor``), see `accessor_catalogue` for its format
and `lb_accessors` for the accessors of the lattice Boltzmann method, which
use the same backends.

Replaces the Jinja templates ``EK_FieldAccessors.tmpl.{h,cuh,cu}``. There are
three namespaces, which exist on both targets:

* ``Scalar``: get, set, add, initialize (set all cells) and add_to_all for a
  scalar field (the density of a species, the electrostatic potential), for a
  cell and a cell interval;
* ``Vector``: the same for a field of a vector with 3 components;
* ``Flux``: get, get_vector and initialize for the flux field, which has 13
  components for each cell (see `flux_layout`). ``get`` returns the 13
  components of a cell, for a cell interval the 13 components of each cell
  in the order of the cells (the old CPU and CUDA accessors returned only
  the first three components for a cell interval: this is fixed).
  ``get_vector`` returns the flux vector of a cell: half of the sum of the
  fluxes through its 26 faces, times the direction of the face, which reads
  the flux field of the neighbors.

The functions of each target and their parameter names are frozen by
``tests/data/ek_field_accessors{,_cuda}_api_*.json``, which were extracted
from the headers that the old template generated.
"""

import pystencils as ps
import sympy as sp

import accessor_backend_cpu
import accessor_backend_cuda
import flux_layout
from accessor_catalogue import (
    Accessor, Catalogue, FieldArg, FieldSpec, Kernel, Kind, Namespace,
    SCALAR, STD_VECTOR, VECTOR3, Value, Variant, array)

_DIM = flux_layout.DIM
_FLUX_COUNT = flux_layout.FLUX_COUNT
_SCALAR = FieldArg("scalar", True, "scalar_field")
_SCALAR_MUT = FieldArg("scalar", False, "scalar_field")
_VEC = FieldArg("vec", True, "vec_field")
_VEC_MUT = FieldArg("vec", False, "vec_field")
_FLUX = FieldArg("flux", True, "flux_field")
_FLUX_MUT = FieldArg("flux", False, "flux_field")


def _targets(auto=False, **cuda_names):
    """Both targets; the CPU returns ``auto`` if requested, the CUDA
    parameters have the names ``cuda_names`` (default name: new name)."""
    return {"cpu": Variant(returns="auto" if auto else None),
            "cuda": Variant(names=cuda_names)}


class EkAccessors:
    """
    Builds the `Catalogue` of the accessors of the electrokinetics with the
    working precision ``dtype`` (a pystencils type).
    """

    def __init__(self, dtype):
        self.dtype = dtype

        # pystencils fields; sizes and strides are set at runtime
        def field(name, rank, fsize):
            return ps.Field.create_generic(
                name, rank, dtype, index_dimensions=1,
                layout="fzyx" if rank == 3 else "zyxf", index_shape=(fsize,))

        self.f_scalar = field("scalar", 3, 1)
        self.f_vec = field("vec", 3, _DIM)
        self.f_flux = field("flux", 3, _FLUX_COUNT)
        self.f_buf = {n: field(f"buf{n}", 3, n)
                      for n in (1, _DIM, _FLUX_COUNT)}
        self.staggered_flux = flux_layout.staggered_flux_field(dtype, "flux")

    # ---- helpers

    def sym(self, name):
        return ps.TypedSymbol(name, self.dtype)

    def buffer(self, n):
        return self.f_buf[n]

    def component_symbols(self, prefix, n):
        """Scalar kernel parameters ``prefix_0`` ... ``prefix_{n-1}``."""
        return [self.sym(f"{prefix}_{i}") for i in range(n)]

    def field_namespace(self, name, field, field_arg, field_arg_mut, n,
                        container, prefix, cuda_cell_names,
                        cuda_interval_names):
        """
        Namespace of a field of ``n`` components of the container
        ``container`` per cell, whose values are called ``prefix``: the
        kernels get, set, add and the two that set the same value
        everywhere, and the accessors for a cell and a cell interval.
        """
        buf = self.buffer(n)
        components = range(n)
        value = Value(container, n, prefix)
        return Namespace(name, kernels=(
            Kernel("get", lambda: [
                ps.Assignment(buf(j), field(j)) for j in components]),
            Kernel("set", lambda: [
                ps.Assignment(field(j), buf(j)) for j in components]),
            Kernel("add", lambda: [
                ps.Assignment(field(j), field(j) + buf(j))
                for j in components]),
            Kernel("broadcast", lambda: [
                ps.Assignment(field(j), s) for j, s in enumerate(
                    self.component_symbols(prefix, n))],
                   components=(prefix, n)),
            Kernel("broadcast_add", lambda: [
                ps.Assignment(field(j), field(j) + s) for j, s in enumerate(
                    self.component_symbols(prefix, n))],
                   components=(prefix, n))),
            accessors=(
            Accessor("get", Kind.CELL, "get", (field_arg,), out=value,
                     targets=_targets(True, **cuda_cell_names)),
            Accessor("set", Kind.CELL, "set", (field_arg_mut,), inp=value,
                     targets=_targets(False, **cuda_cell_names)),
            Accessor("add", Kind.CELL, "add", (field_arg_mut,), inp=value,
                     targets=_targets(False, **cuda_cell_names)),
            Accessor("initialize", Kind.ALL, "broadcast", (field_arg_mut,),
                     inp=value,
                     targets=_targets(False, **cuda_cell_names)),
            Accessor("add_to_all", Kind.ALL, "broadcast_add",
                     (field_arg_mut,), inp=value,
                     targets=_targets(False, **cuda_cell_names)),
            Accessor("get", Kind.INTERVAL, "get", (field_arg,),
                     out=Value(STD_VECTOR, n, "out"),
                     targets=_targets(True, **cuda_interval_names)),
            Accessor("set", Kind.INTERVAL, "set", (field_arg_mut,),
                     inp=Value(STD_VECTOR, n, "values"),
                     targets=_targets(False, **cuda_interval_names))))

    # ---- namespaces

    def scalar(self):
        return self.field_namespace(
            "Scalar", self.f_scalar, _SCALAR, _SCALAR_MUT, 1, SCALAR,
            "value", {"scalar_field": "field"}, {"scalar_field": "field"})

    def vector(self):
        return self.field_namespace(
            "Vector", self.f_vec, _VEC, _VEC_MUT, _DIM, VECTOR3, "vec",
            {"vec_field": "field"}, {})

    def flux(self):
        flux = self.f_flux
        buf_flux, buf_vector = self.buffer(_FLUX_COUNT), self.buffer(_DIM)
        components = range(_FLUX_COUNT)

        def get_vector():
            # half of the sum over the 26 faces of the flux through the face
            # times its direction: the two faces of each pair of opposite
            # directions are added first, so that the pairs of equal and
            # opposite fluxes cancel exactly, in any precision
            pairs = [sp.Symbol(f"pair_{k}") for k in range(_FLUX_COUNT)]
            assignments = []
            total = sp.Matrix([0, 0, 0])
            for pair, direction in zip(pairs, flux_layout.STORED_DIRECTIONS):
                opposite = tuple(-d for d in direction)
                assignments.append(ps.Assignment(
                    pair, flux_layout.lower(
                        self.staggered_flux.face[0, 0, 0, direction], flux)
                    - flux_layout.lower(
                        self.staggered_flux.face[0, 0, 0, opposite], flux)))
                total += pair * sp.Matrix(direction)
            return [*assignments, *(
                ps.Assignment(buf_vector(j), sp.Rational(1, 2) * total[j])
                for j in range(_DIM))]

        values = Value(array(_FLUX_COUNT), _FLUX_COUNT, "values")
        return Namespace("Flux", kernels=(
            Kernel("get", lambda: [
                ps.Assignment(buf_flux(j), flux(j)) for j in components]),
            Kernel("broadcast", lambda: [
                ps.Assignment(flux(j), s) for j, s in enumerate(
                    self.component_symbols("flux", _FLUX_COUNT))],
                   components=("flux", _FLUX_COUNT)),
            Kernel("get_vector", get_vector)),
            accessors=(
            Accessor("get", Kind.CELL, "get", (_FLUX,),
                     out=Value(array(_FLUX_COUNT), _FLUX_COUNT, "value"),
                     targets=_targets(True)),
            Accessor("initialize", Kind.ALL, "broadcast",
                     (_FLUX_MUT,), inp=values,
                     targets=_targets(False, values="flux")),
            Accessor("get", Kind.INTERVAL, "get", (_FLUX,),
                     out=Value(STD_VECTOR, _FLUX_COUNT, "out"),
                     targets=_targets(True)),
            Accessor("get_vector", Kind.CELL, "get_vector", (_FLUX,),
                     out=Value(VECTOR3, _DIM, "result"),
                     targets=_targets(True)),
            Accessor("get_vector", Kind.INTERVAL, "get_vector", (_FLUX,),
                     out=Value(STD_VECTOR, _DIM, "out"),
                     targets=_targets(True))))

    # ---- the catalogue

    def catalogue(self):
        fields = {"scalar": FieldSpec(self.f_scalar, 1),
                  "vec": FieldSpec(self.f_vec, _DIM),
                  "flux": FieldSpec(self.f_flux, _FLUX_COUNT)}
        return Catalogue(
            dtype=self.dtype,
            namespace="walberla::ek::accessor",
            prelude="\nElectrokinetics field accessors.",
            includes={"cpu": ()},
            fields=fields, dim=_DIM,
            namespaces=(self.scalar(), self.vector(), self.flux()),
            component_symbols=self.component_symbols,
            buffer=self.buffer)


def generate_field_accessors(sfg, ctx):
    """
    Add the complete set of field accessors of the electrokinetics to
    ``sfg``: the header-only file ``EKFieldAccessors<Precision>Precision.h``
    on the CPU (``sfg`` must be header-only, see ``KernelFileGenerator``), the
    header ``.cuh`` and the source file ``.cu`` on CUDA.
    """
    catalogue = EkAccessors(ctx.data_type).catalogue()
    backend = (accessor_backend_cuda.CudaBackend if ctx.is_gpu
               else accessor_backend_cpu.CpuBackend)
    backend(sfg, ctx, catalogue).generate()
