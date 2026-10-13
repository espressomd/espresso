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
Catalogue of the LB field accessors of the waLBerla bridge
(``walberla::lbm::accessor``), see `accessor_catalogue` for its format.

Replaces the Jinja templates ``FieldAccessors.tmpl.{h,cuh,cu}``, which were
adapted from the waLBerla ``LatticeModel.tmpl.h``. The namespaces
``Population, Vector, Interpolation, Density, Velocity, Force,
MomentumDensity, PressureTensor`` exist on both targets (the CPU has no
``Interpolation`` and ``Vector::set_from_list``), the CPU has in addition
``EquilibriumDistribution`` and ``Equilibrium``. The functions of each
target and their parameter names are frozen by
``tests/data/field_accessors{,_cuda}_api_*.json``.

The kernels are written as pystencils assignments over fields of the lattice
(``pdfs``, ``velocity``, ``force``, ``vec``) and a buffer through which the
values of the cells are exchanged with the caller. The equations of the LB
method come from `accessor_equations`: an equation set is instantiated with
the reads of the fields (`LbAccessors.reads`). The order of the kernels and
the accessors is the order of the generated CUDA code.
"""

import itertools

import pystencils as ps
import sympy as sp
from pystencils.sympyextensions.typed_sympy import TypeCast

import accessor_backend_cpu
import accessor_backend_cuda
import accessor_equations
from accessor_catalogue import (
    Accessor, Catalogue, CPU_ONLY, CUDA_ONLY, FieldArg, FieldSpec,
    Kernel, Kind, MATRIX3, Namespace, Point, Reduction, ReductionKernel,
    SCALAR, STD_VECTOR, VECTOR3, Value, Variant, array)

# the accessors are written for a three-dimensional lattice
_DIM = 3
_PDF = FieldArg("pdf", True, "pdf_field")
_PDF_MUT = FieldArg("pdf", False, "pdf_field")
_VEL_MUT = FieldArg("vel", False, "velocity_field")
_FORCE = FieldArg("force", True, "force_field")
_FORCE_MUT = FieldArg("force", False, "force_field")
_VEC = FieldArg("vec", True, "vec_field")
_VEC_MUT = FieldArg("vec", False, "vec_field")


def _targets(auto=False, **cuda_names):
    """Both targets; the CPU returns ``auto`` if requested, the CUDA
    parameters have the names ``cuda_names`` (default name: new name)."""
    return {"cpu": Variant(returns="auto" if auto else None),
            "cuda": Variant(names=cuda_names)}


def _instantiate(eqs, inputs, suffix=""):
    """
    Pystencils assignments of the equation set ``eqs``: the input symbols
    (by name) are replaced by ``inputs``, the assigned symbols get ``suffix``
    appended to their names.
    """
    subs = {sp.Symbol(name): value for name, value in inputs.items()}
    assignments = []
    for lhs, rhs in eqs.assignments:
        new = sp.Symbol(lhs.name + suffix)
        assignments.append(ps.Assignment(new, rhs.xreplace(subs)))
        subs[lhs] = new
    return assignments, subs


class LbAccessors:
    """
    Builds the `Catalogue` of the accessors of an lbmpy method with the
    working precision ``dtype`` (a pystencils type).
    """

    def __init__(self, dtype, lb_method):
        self.dtype = dtype
        self.equations = accessor_equations.AccessorEquations(lb_method)
        self.Q = self.equations.Q
        self.D = self.equations.D
        if self.D != _DIM:
            raise ValueError("the field accessors require a "
                             "three-dimensional lattice")
        self.stencil_name = lb_method.stencil.name
        if not self.stencil_name:
            raise ValueError("lb_method uses a stencil that is not supported "
                             "in waLBerla")
        Q, D = self.Q, self.D

        # pystencils fields; sizes and strides are set at runtime
        def field(name, rank, fsize, dtype=None, **kwargs):
            return ps.Field.create_generic(
                name, rank, dtype or self.dtype, index_dimensions=1,
                layout="fzyx" if rank == 3 else "zyxf",
                index_shape=(fsize,), **kwargs)

        custom = dict(field_type=ps.FieldType.CUSTOM)
        self.f_pdf = field("pdfs", 3, Q)
        self.f_vel = field("velocity", 3, D)
        self.f_force = field("force", 3, D)
        self.f_vec = field("vec", 3, D)
        # fields of the kernels over lists, which address the lattice fields
        # with absolute indices
        self.f_pdf_abs = field("pdfs", 3, Q, **custom)
        self.f_vec_abs = field("vec", 3, D, **custom)
        self.f_pos = field("pos", 1, D)
        self.f_indices = field("indices", 1, D, dtype=ps.create_type("int32"))
        self.f_values = field("values", 1, D)
        self.f_out = {n: field(f"out{n}", 1, n) for n in (1, D)}
        self.f_buf = {n: field(f"buf{n}", 3, n) for n in (1, D, D * D, Q)}

    # ---- helpers

    def sym(self, name):
        return ps.TypedSymbol(name, self.dtype)

    def buffer(self, n):
        return self.f_buf[n]

    def component_symbols(self, prefix, n):
        """Scalar kernel parameters ``prefix_0`` ... ``prefix_{n-1}``."""
        return [self.sym(f"{prefix}_{i}") for i in range(n)]

    def reads(self, eqs, pdf=None, force=None, **others):
        """Replacement of the input symbols of an equation set: the
        populations ``f_i`` by ``pdf(i)``, the forces by ``force(i)``."""
        E = self.equations
        inputs = dict(others)
        if pdf is not None:
            inputs.update({E.pdf_names[i]: pdf(i) for i in E.used_pdfs(eqs)})
        if force is not None:
            inputs.update({E.force_names[i]: force(i)
                           for i in E.used_forces(eqs)})
        return inputs

    def equilibrium(self, pdf, rho, u):
        """Assignments that set all populations of ``pdf`` to the
        equilibrium of the density ``rho`` and the velocity ``u``."""
        E = self.equations
        delta_rho = sp.Symbol("delta_rho_eq")
        subs = {sp.Symbol("rho"): rho, sp.Symbol("delta_rho"): delta_rho,
                **{sp.Symbol(n): v for n, v in zip(E.u_names, u)}}
        return [ps.Assignment(delta_rho, rho - 1)] + [
            ps.Assignment(pdf(i), expr.xreplace(subs))
            for i, expr in enumerate(E.equilibrium)]

    def momenta_over_density(self, rho_inv):
        """Velocity components ``md_j * rho_inv``."""
        return [sp.Symbol(f"md_{j}") * rho_inv for j in range(self.D)]

    def interpolation_weights(self):
        """
        Assignments of the trilinear interpolation weights of the particle
        position ``pos``, which is relative to the origin of the local
        halo-free domain, with cells centered at half-integer positions.

        The same weights are computed by the C++ function
        ``calculate_weights`` of the CUDA backend (``Interpolation::
        add_force``): keep both in sync (they are covered by the
        interpolation tests and the GPU equivalence harness).
        Returns the assignments, the corner cell (in ``pdf`` index space,
        which starts ``gl`` cells below the origin), and a function that
        returns the weight of the corner ``(i, j, k)``.
        """
        half = sp.Rational(1, 2)
        int64 = ps.create_type("int64")
        gl = ps.TypedSymbol("ghost_layers", int64)
        assignments, corner, weights = [], [], []
        for d in range(self.D):
            fractional = sp.Symbol(f"fractional_index_{d}")
            nmp = sp.Symbol(f"nmp_{d}")
            distance = sp.Symbol(f"distance_{d}")
            c = ps.TypedSymbol(f"corner_{d}", int64)
            assignments += [
                ps.Assignment(fractional, self.f_pos(d) - half),
                ps.Assignment(nmp, sp.floor(fractional)),
                ps.Assignment(distance, fractional - nmp - half),
                ps.Assignment(c, TypeCast(nmp, int64) + gl)]
            corner.append(c)
            weights.append((half - distance, half + distance))

        def weight(i, j, k):
            return weights[0][i] * weights[1][j] * weights[2][k]
        return assignments, corner, weight

    def interpolate(self, fsize, out, cell_values):
        """
        Assignments that interpolate ``fsize`` values per particle,
        returned by ``cell_values(corner, n)`` for the corner cell with
        absolute indices ``corner`` and the unique name suffix ``n``.
        The corners are accumulated in the order of the nested loops
        ``i, j, k``.
        """
        assignments, corner, weight = self.interpolation_weights()
        accumulated = [None] * fsize
        for n, (i, j, k) in enumerate(itertools.product((0, 1), repeat=3)):
            cell = (corner[0] + i, corner[1] + j, corner[2] + k)
            w = sp.Symbol(f"weight_{n}")
            assignments.append(ps.Assignment(w, weight(i, j, k)))
            prologue, values = cell_values(cell, n)
            assignments += prologue
            for cf in range(fsize):
                total = w * values[cf]
                if accumulated[cf] is not None:
                    total = accumulated[cf] + total
                accumulated[cf] = sp.Symbol(f"accumulated_{cf}_{n}")
                assignments.append(ps.Assignment(accumulated[cf], total))
        assignments += [ps.Assignment(out(cf), accumulated[cf])
                        for cf in range(fsize)]
        return assignments

    # ---- the catalogue

    def catalogue(self):
        fields = {
            "pdf": FieldSpec(self.f_pdf, self.Q, self.f_pdf_abs),
            "vel": FieldSpec(self.f_vel, self.D),
            "force": FieldSpec(self.f_force, self.D),
            "vec": FieldSpec(self.f_vec, self.D, self.f_vec_abs),
            "pos": FieldSpec(self.f_pos, self.D),
            "indices": FieldSpec(self.f_indices, self.D),
            "values": FieldSpec(self.f_values, self.D),
            **{f"out{n}": FieldSpec(f, n) for n, f in self.f_out.items()}}
        return Catalogue(
            dtype=self.dtype,
            namespace="walberla::lbm::accessor",
            prelude="\nLattice field accessors.\n"
                    "Adapted from the waLBerla source file "
                    "python/lbmpy_walberla/templates/LatticeModel.tmpl.h "
                    "(commit a16141524c58ab88386e2a0f8fdd7c63c5edd704).",
            includes={"cpu": (f"stencil/{self.stencil_name}.h",)},
            fields=fields, dim=self.D,
            namespaces=(
                self.population(), self.vector(),
                self.equilibrium_distribution(), self.equilibrium_ns(),
                self.interpolation(), self.density(), self.velocity(),
                self.force(), self.momentum_density(),
                self.pressure_tensor()),
            component_symbols=self.component_symbols,
            buffer=self.buffer)

    # ---- namespace Population

    def population(self):
        Q = self.Q
        E = self.equations
        f_pdf, f_vel, f_force, buf = (self.f_pdf, self.f_vel, self.f_force,
                                      self.buffer(Q))
        md = E.momentum_density
        rho, rho_inv = sp.Symbol("rho"), sp.Symbol("rho_inv")

        def get():
            return [ps.Assignment(buf(i), f_pdf(i)) for i in range(Q)]

        def set_():
            return [ps.Assignment(f_pdf(i), buf(i)) for i in range(Q)]

        def broadcast():
            return [ps.Assignment(f_pdf(i), s) for i, s in enumerate(
                self.component_symbols("pop", Q))]

        def set_vel():
            # the populations are read from the buffer
            temporaries, _ = _instantiate(
                md, self.reads(md, pdf=buf, force=f_force.center))
            return [
                *temporaries, ps.Assignment(rho_inv, 1 / rho),
                *(ps.Assignment(f_pdf(i), buf(i)) for i in range(Q)),
                *(ps.Assignment(f_vel(j), m) for j, m in enumerate(
                    self.momenta_over_density(rho_inv)))]

        pop = Value(array(Q), Q, "pop")
        values = Value(STD_VECTOR, Q, "values")
        out = Value(STD_VECTOR, Q, "out")
        pdf_vel_force = (_PDF_MUT, _VEL_MUT, _FORCE)
        return Namespace("Population", kernels=(
            Kernel("get", get), Kernel("set", set_),
            Kernel("broadcast", broadcast, components=("pop", Q)),
            Kernel("set_vel", set_vel)), accessors=(
            Accessor("get", Kind.CELL, "get", (_PDF,),
                     out=Value(array(Q), Q, "pop"),
                     targets=_targets(auto=True)),
            Accessor("set", Kind.CELL, "set", (_PDF_MUT,), inp=pop),
            Accessor("set", Kind.CELL, "set_vel", pdf_vel_force, inp=pop),
            Accessor("initialize", Kind.ALL, "broadcast", (_PDF_MUT,),
                     inp=pop),
            Accessor("get", Kind.INTERVAL, "get", (_PDF,), out=out,
                     targets=_targets(auto=True)),
            Accessor("set", Kind.INTERVAL, "set", (_PDF_MUT,), inp=values),
            Accessor("set", Kind.INTERVAL, "set_vel", pdf_vel_force,
                     inp=values)))

    # ---- namespace Vector

    def vector(self):
        D = self.D
        f_vec, buf = self.f_vec, self.buffer(D)
        comps = range(D)

        def get():
            return [ps.Assignment(buf(j), f_vec(j)) for j in comps]

        def set_():
            return [ps.Assignment(f_vec(j), buf(j)) for j in comps]

        def add():
            return [ps.Assignment(f_vec(j), f_vec(j) + buf(j))
                    for j in comps]

        def broadcast():
            return [ps.Assignment(f_vec(j), s) for j, s in enumerate(
                self.component_symbols("vec", D))]

        def broadcast_add():
            return [ps.Assignment(f_vec(j), f_vec(j) + s)
                    for j, s in enumerate(self.component_symbols("vec", D))]

        def set_from_list():
            idx = tuple(self.f_indices(j) for j in comps)
            return [ps.Assignment(
                self.f_vec_abs.absolute_access(idx, (j,)), self.f_values(j))
                for j in comps]

        vec = Value(VECTOR3, D, "vec")
        # the names of the parameters of the CUDA functions
        cell_names = _targets(auto=True, vec_field="field")
        cell_names_void = _targets(vec_field="field")
        return Namespace("Vector", kernels=(
            Kernel("get", get), Kernel("set", set_), Kernel("add", add),
            Kernel("broadcast", broadcast, components=("vec", D)),
            Kernel("broadcast_add", broadcast_add,
                   components=("vec", D))), accessors=(
            Accessor("get", Kind.CELL, "get", (_VEC,), out=vec,
                     targets=cell_names),
            Accessor("set", Kind.CELL, "set", (_VEC_MUT,), inp=vec,
                     targets=cell_names_void),
            Accessor("add", Kind.CELL, "add", (_VEC_MUT,), inp=vec,
                     targets=cell_names_void),
            Accessor("initialize", Kind.ALL, "broadcast", (_VEC_MUT,),
                     inp=vec, targets=cell_names_void),
            Accessor("add_to_all", Kind.ALL, "broadcast_add", (_VEC_MUT,),
                     inp=vec, targets=cell_names_void),
            Accessor("get", Kind.INTERVAL, "get", (_VEC,),
                     out=Value(STD_VECTOR, D, "out"),
                     targets=_targets(auto=True)),
            Accessor("set", Kind.INTERVAL, "set", (_VEC_MUT,),
                     inp=Value(STD_VECTOR, D, "values")),
            Accessor("set_from_list", Kind.LIST,
                     Kernel("set_from_list", set_from_list),
                     (FieldArg("vec", True, "field"),), dims=D,
                     targets=CUDA_ONLY)))

    # ---- namespace EquilibriumDistribution

    def equilibrium_distribution(self):
        E = self.equations
        point = Point(tuple(E.u_names), tuple(E.equilibrium_from_direction))
        return Namespace("EquilibriumDistribution", kernels=(), accessors=(
            Accessor("get", Kind.POINT, point=point, targets=CPU_ONLY),))

    # ---- namespace Equilibrium

    def equilibrium_ns(self):
        D = self.D
        f_pdf, buf = self.f_pdf, self.buffer(self.D)

        def set_():
            # the velocity is read once, not in every equilibrium
            u = [sp.Symbol(f"u_eq_{j}") for j in range(D)]
            return [*(ps.Assignment(uj, buf(j)) for j, uj in enumerate(u)),
                    *self.equilibrium(f_pdf, self.sym("rho"), u)]

        return Namespace("Equilibrium", kernels=(
            Kernel("set", set_, targets=("cpu",)),), accessors=(
            Accessor("set", Kind.CELL, "set", (_PDF_MUT,),
                     inp=Value(VECTOR3, D, "u"), scalars=("rho",),
                     targets=CPU_ONLY),))

    # ---- namespace Interpolation

    def interpolation(self):
        D = self.D
        E = self.equations
        density_sym = self.sym("density")
        scaled = E.density_scaled

        def cell_density(cell, n):
            reads = self.reads(
                scaled, pdf=lambda q: self.f_pdf_abs.absolute_access(
                    cell, (q,)), density=density_sym)
            temporaries, subs = _instantiate(scaled, reads, f"_{n}")
            return temporaries, [subs[sp.Symbol("rho")]]

        def cell_velocity(cell, n):
            return [], [self.f_vec_abs.absolute_access(cell, (cf,))
                        for cf in range(D)]

        def get_rho():
            return self.interpolate(1, self.f_out[1], cell_density)

        def get_vel():
            return self.interpolate(D, self.f_out[D], cell_velocity)

        vec_out = Value(STD_VECTOR, 1, "out")
        return Namespace("Interpolation", kernels=(
            Kernel("get_rho", get_rho, targets=("cuda",)),
            Kernel("get_vel", get_vel, targets=("cuda",))), accessors=(
            Accessor("get_rho", Kind.GATHER, "get_rho",
                     (FieldArg("pdf", True, "field"),), out=vec_out,
                     scalars=("density",), dims=1, targets=CUDA_ONLY),
            Accessor("get_vel", Kind.GATHER, "get_vel",
                     (FieldArg("vec", True, "field"),),
                     out=Value(STD_VECTOR, D, "out"), dims=D,
                     targets=CUDA_ONLY),
            Accessor("add_force", Kind.SCATTER, None,
                     (FieldArg("vec", True, "field"),), dims=D,
                     targets=CUDA_ONLY)))

    # ---- namespace Density

    def density(self):
        D = self.D
        E = self.equations
        f_pdf, buf = self.f_pdf, self.buffer(1)
        density_sym = self.sym("density")
        rho, rho_inv = sp.Symbol("rho"), sp.Symbol("rho_inv")
        scaled = E.density_scaled
        unshifted = E.unshifted_momentum_density

        def get():
            temporaries, _ = _instantiate(
                scaled, self.reads(scaled, pdf=f_pdf, density=density_sym))
            return [*temporaries, ps.Assignment(buf(0), rho)]

        def set_():
            # the density is changed at constant velocity (without the half
            # force shift)
            temporaries, _ = _instantiate(
                unshifted, self.reads(unshifted, pdf=f_pdf))
            rho_new = sp.Symbol("rho_new")
            u_old = [sp.Symbol(f"u_old_{j}") for j in range(D)]
            return [
                *temporaries, ps.Assignment(rho_inv, 1 / rho),
                *(ps.Assignment(u, m) for u, m in zip(
                    u_old, self.momenta_over_density(rho_inv))),
                ps.Assignment(rho_new, buf(0) / density_sym),
                *self.equilibrium(f_pdf, rho_new, u_old)]

        return Namespace("Density", kernels=(
            Kernel("get", get), Kernel("set", set_)), accessors=(
            Accessor("get", Kind.CELL, "get", (_PDF,),
                     out=Value(SCALAR, 1, "rho"), scalars=("density",)),
            Accessor("get", Kind.INTERVAL, "get", (_PDF,),
                     out=Value(STD_VECTOR, 1, "out"), scalars=("density",)),
            Accessor("set", Kind.CELL, "set", (_PDF_MUT,),
                     inp=Value(SCALAR, 1, "rho_in"), scalars=("density",),
                     targets=_targets(rho_in="rho")),
            Accessor("set", Kind.INTERVAL, "set", (_PDF_MUT,),
                     inp=Value(STD_VECTOR, 1, "values"),
                     scalars=("density",))))

    # ---- namespace Velocity

    def velocity(self):
        D = self.D
        E = self.equations
        f_pdf, f_vel, f_force, buf = (self.f_pdf, self.f_vel, self.f_force,
                                      self.buffer(D))
        rho, rho_inv = sp.Symbol("rho"), sp.Symbol("rho_inv")
        md = E.momentum_density
        dens = E.density
        setter = E.velocity_setter

        def get():
            temporaries, _ = _instantiate(
                md, self.reads(md, pdf=f_pdf, force=f_force.center))
            return [*temporaries, ps.Assignment(rho_inv, 1 / rho),
                    *(ps.Assignment(buf(j), m) for j, m in enumerate(
                        self.momenta_over_density(rho_inv)))]

        def set_():
            # the velocity is shifted by half the force, the stored
            # velocity is not
            temporaries, _ = _instantiate(
                dens, self.reads(dens, pdf=f_pdf))
            shift, _ = _instantiate(
                setter, self.reads(setter, force=f_force.center,
                                   **{n: buf(j) for j, n in
                                      enumerate(E.u_in_names)}))
            return [*temporaries, *shift,
                    *(ps.Assignment(f_vel(j), buf(j)) for j in range(D)),
                    *self.equilibrium(
                        f_pdf, rho, [sp.Symbol(n) for n in E.u_names])]

        pdf_vel_force = (_PDF_MUT, _VEL_MUT, _FORCE)
        return Namespace("Velocity", kernels=(
            Kernel("get", get), Kernel("set", set_)), accessors=(
            Accessor("get", Kind.CELL, "get", (_PDF, _FORCE),
                     out=Value(VECTOR3, D, "vec"),
                     targets=_targets(auto=True)),
            # the buffer holds the vectors of all cells of the interval
            Accessor("get", Kind.INTERVAL, "get", (_PDF, _FORCE),
                     out=Value(STD_VECTOR, D, "out"),
                     targets=_targets(auto=True)),
            Accessor("set", Kind.CELL, "set", pdf_vel_force,
                     inp=Value(VECTOR3, D, "u")),
            Accessor("set", Kind.INTERVAL, "set", pdf_vel_force,
                     inp=Value(STD_VECTOR, D, "values"))))

    # ---- namespace Force

    def force(self):
        D = self.D
        E = self.equations
        f_pdf, f_vel, f_force, buf = (self.f_pdf, self.f_vel, self.f_force,
                                      self.buffer(D))
        density_sym = self.sym("density")
        rho, rho_inv = sp.Symbol("rho"), sp.Symbol("rho_inv")
        eqs = E.momentum_density_force_setter

        def set_():
            temporaries, _ = _instantiate(
                eqs, self.reads(
                    eqs, pdf=f_pdf, density=density_sym,
                    **{n: buf(j) for j, n in enumerate(E.force_names)}))
            density_inv = sp.Symbol("density_inv")
            return [
                *temporaries, ps.Assignment(rho_inv, 1 / rho),
                ps.Assignment(density_inv, 1 / density_sym),
                *(ps.Assignment(f_force(j), buf(j) * density_inv)
                  for j in range(D)),
                *(ps.Assignment(f_vel(j), m) for j, m in enumerate(
                    self.momenta_over_density(rho_inv)))]

        fields = (_PDF, _VEL_MUT, _FORCE_MUT)
        return Namespace("Force", kernels=(Kernel("set", set_),), accessors=(
            Accessor("set", Kind.CELL, "set", fields,
                     inp=Value(VECTOR3, D, "force"), scalars=("density",),
                     targets=_targets(force="u")),
            Accessor("set", Kind.INTERVAL, "set", fields,
                     inp=Value(STD_VECTOR, D, "values"),
                     scalars=("density",))))

    # ---- namespace MomentumDensity

    def momentum_density(self):
        D = self.D
        E = self.equations
        f_pdf, f_force = self.f_pdf, self.f_force
        density_sym = self.sym("density")
        # only the shifted momentum density is needed: this skips the
        # density and the populations it depends on
        eqs = E.momentum_density.pruned([f"md_{j}" for j in range(D)])

        def reduce():
            temporaries, _ = _instantiate(
                eqs, self.reads(eqs, pdf=f_pdf, force=f_force.center))
            return temporaries, [sp.Symbol(f"md_{j}") * density_sym
                                 for j in range(D)]

        # the kernel is named ``get`` in the generated CUDA code
        return Namespace("MomentumDensity", kernels=(
            ReductionKernel("get", reduce),), accessors=(
            Accessor("reduce", Kind.REDUCE, "get", (_PDF, _FORCE),
                     scalars=("density",),
                     reduction=Reduction(VECTOR3, "momentum_density"),
                     targets=_targets(auto=True)),))

    # ---- namespace PressureTensor

    def pressure_tensor(self):
        D = self.D
        E = self.equations
        f_pdf, buf = self.f_pdf, self.buffer(D * D)
        density_sym = self.sym("density")
        eqs = E.second_moment
        offsets = E.inverse_neighbor_offsets
        comps = range(D * D)
        p = [sp.Symbol(f"p_{i}") for i in comps]

        def temporaries():
            # the populations are pulled from the neighbor cells
            used = {E.pdf_names[i]: f_pdf[offsets[i]](i)
                    for i in E.used_pdfs(eqs)}
            return _instantiate(eqs, used)[0]

        def get():
            return [*temporaries(), *(ps.Assignment(buf(i), p[i] * density_sym)
                                      for i in comps)]

        def reduce():
            return temporaries(), list(p)

        return Namespace("PressureTensor", kernels=(
            Kernel("get", get), ReductionKernel("reduce", reduce)),
            accessors=(
            Accessor("get", Kind.CELL, "get", (_PDF,),
                     out=Value(MATRIX3, D * D, "out"),
                     scalars=("density",), targets=_targets(auto=True)),
            Accessor("get", Kind.INTERVAL, "get", (_PDF,),
                     out=Value(STD_VECTOR, D * D, "out"),
                     scalars=("density",), targets=_targets(auto=True)),
            Accessor("reduce", Kind.REDUCE, "reduce", (_PDF,),
                     scalars=("density",),
                     reduction=Reduction(MATRIX3, "pressure_tensor",
                                         "density"),
                     targets=_targets(auto=True))))


def generate_field_accessors(sfg, ctx, lb_method):
    """
    Add the complete set of field accessors of ``lb_method`` to ``sfg``:
    the header-only file ``FieldAccessors<Precision>Precision.h`` on the CPU
    (``sfg`` must be header-only, see ``KernelFileGenerator``), the header
    ``.cuh`` and the source file ``.cu`` on CUDA.
    """
    catalogue = LbAccessors(ctx.data_type, lb_method).catalogue()
    backend = (accessor_backend_cuda.CudaBackend if ctx.is_gpu
               else accessor_backend_cpu.CpuBackend)
    backend(sfg, ctx, catalogue).generate()
