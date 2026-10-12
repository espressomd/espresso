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
CPU lattice field accessors of the waLBerla bridge, generated with the
pystencils-sfg composer. Replaces the Jinja template ``FieldAccessors.tmpl.h``
(adapted from the waLBerla ``LatticeModel.tmpl.h``).

The header-only file defines inline functions in the namespaces
``walberla::lbm::accessor::{Population, Vector, EquilibriumDistribution,
Equilibrium, Density, Velocity, Force, MomentumDensity, PressureTensor}``.
Each function is added with ``sfg.function``; the statements of the
function bodies are C++ strings assembled in Python. The equations of the
LB method are printed by `accessor_equations`.
"""

from pystencils.types import PsCustomType, constify
from pystencilssfg.composer.basic_composer import (
    SfgFunctionSequencer, make_sequence)
from pystencilssfg.ir import SfgEntityDef, SfgFunction
from pystencilssfg.lang import AugExpr, cpptype
from pystencilssfg.lang.cpp import std
from sweepgen.api import (
    Cell, CellInterval, Direction, GhostLayerFieldPtr, Vector3)

import accessor_equations

# the accessors are written for a three-dimensional lattice
_DIM = 3

_Array = cpptype("std::array< {T}, {N}u >", "<array>")
_Matrix3 = cpptype("walberla::Matrix3< {T} >", "core/math/Matrix3.h")

AUTO = PsCustomType("auto")


class _OverloadSequencer(SfgFunctionSequencer):
    """
    Function sequencer that allows overloads. ``SfgFunctionSequencer``
    registers every function by name in its namespace and rejects a second
    function of the same name, which the accessors need (e.g. the ``get``
    functions for a cell and for a cell interval). The overloads are
    written to the header without registering them.
    """

    def __call__(self, *args):
        tree = make_sequence(*args)
        func = SfgFunction(
            self._name, self._cursor.current_namespace, tree,
            return_type=self._return_type, inline=True,
            constexpr=self._constexpr, attributes=self._attributes,
            required_params=self._params)
        self._cursor.write_header(SfgEntityDef(func))


def _ui(i):
    """Spelling of an index literal: ``uint_t{3u}``."""
    return f"uint_t{{{i}u}}"


def _interval_loops(*body):
    """
    Triple loop over the cells of the ``CellInterval ci``, as arguments of a
    function sequence (loop headers followed by blocks).
    """
    return ("for (auto x = ci.xMin(); x <= ci.xMax(); ++x)",
            ("for (auto y = ci.yMin(); y <= ci.yMax(); ++y)",
             ("for (auto z = ci.zMin(); z <= ci.zMax(); ++z)", body)))


def _interior_loops(*body):
    """Triple loop over the interior of ``pdf_field`` (z outermost)."""
    return ("for (cell_idx_t z = 0; z < cell_idx_c(pdf_field->zSize()); ++z)",
            ("for (cell_idx_t y = 0; y < cell_idx_c(pdf_field->ySize()); ++y)",
             ("for (cell_idx_t x = 0; x < cell_idx_c(pdf_field->xSize()); ++x)", body)))


def _all_cells_macro(field, *body):
    """Loop over all cells including ghost layers (waLBerla macro)."""
    return (f"WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ({field}, {{",
            *body, "});")


class _FieldAccessors:
    def __init__(self, sfg, ctx, lb_method):
        self.sfg = sfg
        self.dtype = ctx.data_type
        self.T = self.dtype.c_string().strip()
        self.equations = accessor_equations.AccessorEquations(
            lb_method, self.dtype)
        self.Q = self.equations.Q
        self.D = self.equations.D
        assert self.D == _DIM
        self.stencil_name = lb_method.stencil.name
        if not self.stencil_name:
            raise ValueError("lb_method uses a stencil that is not supported "
                             "in waLBerla")

        # parameter types, as in the legacy header
        T, Q, D = self.dtype, self.Q, self.D
        self.pdf_const = GhostLayerFieldPtr(
            T, D, Q, const=True).var("pdf_field")
        self.pdf_mut = GhostLayerFieldPtr(T, D, Q).var("pdf_field")
        # an AugExpr can be bound to a name only once: these are factories
        self.vec_const = lambda name: GhostLayerFieldPtr(
            T, D, D, const=True).var(name)
        self.vec_mut = lambda name: GhostLayerFieldPtr(T, D, D).var(name)
        self.cell = Cell(const=True, ref=True).var("cell")
        self.ci = CellInterval(const=True, ref=True).var("ci")
        self.vector = lambda name: Vector3(T, const=True, ref=True).var(name)
        self.density_param = sfg.var("density", constify(T))
        self.rho_param = sfg.var("rho", constify(T))
        self.values = std.vector(T, const=True, ref=True).var("values")
        self.pop = AugExpr(_Array(T=T, N=Q, const=True, ref=True)).var("pop")
        self.vector_return = std.vector(T).dtype
        self.vector3_return = Vector3(T).dtype
        self.matrix3_return = _Matrix3(T=T)
        self.array_return = _Array(T=T, N=Q)

    # ---- statements shared by several accessors

    def pdf_loads(self, indices, field="pdf_field", ptr="xyz0"):
        return [f"const {self.T} f_{i} = {field}->getF(&{ptr}, {_ui(i)});"
                for i in indices]

    def force_loads(self, eqs):
        """Bind the force symbols of ``eqs`` to the cell ``(x, y, z)``."""
        return [f"const {self.T} force_{i} = "
                f"force_field->get(x, y, z, {_ui(i)});"
                for i in self.equations.used_forces(eqs)]

    def bind_vector(self, names, vector, eqs=None):
        """Bind scalar symbols to the components of a vector or pointer."""
        return [f"const {self.T} {name} = {vector}[{_ui(i)}];"
                for i, name in enumerate(names)
                if eqs is None or name in eqs.inputs]

    def cell_coordinates(self):
        return [f"const auto {c} = cell.{c}();" for c in "xyz"]

    def pdf_offset_loads(self, indices):
        """Load populations pulled from the neighbor cells ``-c_i``."""
        lines = []
        offsets = self.equations.inverse_neighbor_offsets
        for i in indices:
            coords = []
            for c, o in zip("xyz", offsets[i]):
                coords.append(c if o == 0 else f"{c} {'+' if o > 0 else '-'} "
                              f"{abs(o)}")
            lines.append(f"const {self.T} f_{i} = pdf_field->get("
                         f"{', '.join(coords)}, {_ui(i)});")
        return lines

    def equilibrium_call(self, cell_expr):
        u = ", ".join(self.equations.u_names)
        return (f"Equilibrium::set(pdf_field, Vector{self.D}<{self.T}>({u}), "
                f"rho, {cell_expr});")

    def fn(self, name, *params, returns=None):
        """Start an inline function; the sequencer is called with the body."""
        sequencer = _OverloadSequencer(self.sfg._cursor, name)
        if returns is not None:
            sequencer = sequencer.returns(returns)
        return sequencer.params(*params)

    # ---- generate

    def generate(self):
        sfg = self.sfg
        for header in ("core/DataTypes.h", "core/Abort.h",
                       "core/math/Matrix3.h",
                       f"stencil/{self.stencil_name}.h",
                       "field/iterators/IteratorMacros.h",
                       "<array>", "<cassert>", "<iterator>", "<vector>"):
            sfg.include(header)
        sfg.prelude("\nLattice field accessors.\n"
                    "Adapted from the waLBerla source file "
                    "python/lbmpy_walberla/templates/LatticeModel.tmpl.h "
                    "(commit a16141524c58ab88386e2a0f8fdd7c63c5edd704).")
        sfg.namespace("walberla::lbm::accessor")
        for name in ("population", "vector", "equilibrium_distribution",
                     "equilibrium", "density", "velocity", "force",
                     "momentum_density", "pressure_tensor"):
            getattr(self, f"gen_{name}")()

    # ---- namespace Population

    def gen_population(self):
        T, Q, D = self.T, self.Q, self.D
        E = self.equations
        pdf, pdf_mut, cell, ci = (self.pdf_const, self.pdf_mut, self.cell,
                                  self.ci)
        pop, values = self.pop, self.values
        vel = self.vec_mut("velocity_field")
        frc = self.vec_const("force_field")
        all_q = range(Q)
        with self.sfg.namespace("Population"):
            self.fn("get", pdf, cell, returns=AUTO)(
                f"{T} const &xyz0 = pdf_field->get(cell, {_ui(0)});",
                f"std::array<{T}, {Q}u> pop;",
                *(f"pop[{i}u] = pdf_field->getF(&xyz0, {_ui(i)});"
                  for i in all_q),
                "return pop;")

            self.fn("set", pdf_mut, pop, cell)(
                f"{T} &xyz0 = pdf_field->get(cell, {_ui(0)});",
                *(f"pdf_field->getF(&xyz0, {_ui(i)}) = pop[{i}u];"
                  for i in all_q))

            def set_with_velocity(eqs, pop_expr):
                used = set(E.used_pdfs(eqs))
                lines = []
                for i in all_q:
                    store = (f"pdf_field->getF(&xyz0, {_ui(i)}) = "
                             f"{pop_expr}[{i}u];")
                    lines.append(f"const {T} f_{i} = {store}" if i in used
                                 else store)
                return lines

            eqs = E.momentum_density
            self.fn("set", pdf_mut, vel, frc, pop, cell)(
                f"auto &xyz0 = pdf_field->get(cell, {_ui(0)});",
                *set_with_velocity(eqs, "pop"),
                *self.cell_coordinates(),
                *self.force_loads(eqs),
                *eqs.lines(),
                f"const auto rho_inv = {T}{{1}} / rho;",
                *(f"velocity_field->get(cell, {_ui(i)}) = md_{i} * rho_inv;"
                  for i in range(D)))

            self.fn("initialize", pdf_mut,
                    AugExpr(_Array(T=self.dtype, N=Q, const=True, ref=True)
                            ).var("pop"))(
                *_all_cells_macro(
                    "pdf_field",
                    f"{T} &xyz0 = pdf_field->get(x, y, z, {_ui(0)});",
                    *(f"pdf_field->getF(&xyz0, {_ui(i)}) = pop[{i}u];"
                      for i in all_q)))

            self.fn("get", pdf, ci, returns=AUTO)(
                f"std::vector<{T}> out;",
                f"out.reserve(ci.numCells() * uint_t({Q}u));",
                *_interval_loops(
                    f"{T} const &xyz0 = pdf_field->get(x, y, z, {_ui(0)});",
                    *(f"out.emplace_back(pdf_field->getF(&xyz0, {_ui(i)}));"
                      for i in all_q)),
                "return out;")

            self.fn("set", pdf_mut, values, ci)(
                f"assert(uint_c(values.size()) == ci.numCells() * "
                f"uint_t({Q}u));",
                "auto pop = values.data();",
                *_interval_loops(
                    f"{T} &xyz0 = pdf_field->get(x, y, z, {_ui(0)});",
                    *(f"pdf_field->getF(&xyz0, {_ui(i)}) = pop[{i}u];"
                      for i in all_q),
                    f"std::advance(pop, {Q});"))

            self.fn("set", pdf_mut, vel, frc, values, ci)(
                f"assert(uint_c(values.size()) == ci.numCells() * "
                f"uint_t({Q}u));",
                "auto pop = values.data();",
                *_interval_loops(
                    f"{T} &xyz0 = pdf_field->get(x, y, z, {_ui(0)});",
                    *set_with_velocity(eqs, "pop"),
                    *self.force_loads(eqs),
                    *eqs.lines(),
                    f"const auto rho_inv = {T}{{1}} / rho;",
                    *(f"velocity_field->get(x, y, z, {_ui(i)}) = "
                      f"md_{i} * rho_inv;" for i in range(D)),
                    f"std::advance(pop, {Q});"))

    # ---- namespace Vector

    def gen_vector(self):
        T, D = self.T, self.D
        vec_c, vec_m = self.vec_const("vec_field"), \
            self.vec_mut("vec_field")
        vec = self.vector("vec")
        cell, ci, values = self.cell, self.ci, self.values
        comps = range(D)
        with self.sfg.namespace("Vector"):
            self.fn("get", vec_c, cell, returns=AUTO)(
                f"const {T} &xyz0 = vec_field->get(cell, {_ui(0)});",
                f"Vector{D}<{T}> vec;",
                *(f"vec[{i}] = vec_field->getF(&xyz0, {_ui(i)});"
                  for i in comps),
                "return vec;")

            for name, op in (("set", "="), ("add", "+=")):
                self.fn(name, vec_m, vec, cell)(
                    f"{T} &xyz0 = vec_field->get(cell, {_ui(0)});",
                    *(f"vec_field->getF(&xyz0, {_ui(i)}) {op} vec[{i}u];"
                      for i in comps))

            for name, op in (("initialize", "="), ("add_to_all", "+=")):
                self.fn(name, vec_m, vec)(
                    *_all_cells_macro(
                        "vec_field",
                        f"{T} &xyz0 = vec_field->get(x, y, z, {_ui(0)});",
                        *(f"vec_field->getF(&xyz0, {_ui(i)}) {op} "
                          f"vec[{i}u];" for i in comps)))

            self.fn("get", vec_c, ci, returns=AUTO)(
                f"std::vector<{T}> out;",
                f"out.reserve(ci.numCells() * uint_t({D}u));",
                *_interval_loops(
                    f"const {T} &xyz0 = vec_field->get(x, y, z, {_ui(0)});",
                    *(f"out.emplace_back(vec_field->getF(&xyz0, {_ui(i)}));"
                      for i in comps)),
                "return out;")

            self.fn("set", vec_m, values, ci)(
                f"assert(uint_c(values.size()) == ci.numCells() * "
                f"uint_t({D}u));",
                "auto values_ptr = values.data();",
                *_interval_loops(
                    f"{T} &xyz0 = vec_field->get(x, y, z, {_ui(0)});",
                    *(f"vec_field->getF(&xyz0, {_ui(i)}) = values_ptr[{i}u];"
                      for i in comps),
                    f"std::advance(values_ptr, {D});"))

    # ---- namespace EquilibriumDistribution

    def gen_equilibrium_distribution(self):
        T, D = self.T, self.D
        E = self.equations
        direction = Direction(const=True).var("direction")
        u = self.vector("u")
        rho = AugExpr(self.dtype).var("rho")
        switch = self.sfg.switch("direction", autobreak=False)
        for name, expr in E.equilibrium_from_direction:
            switch.case(f"stencil::{name}")(
                f"return {accessor_equations.print_expression(expr, self.dtype)};")
        switch.default('WALBERLA_ABORT("Invalid Direction")')
        with self.sfg.namespace("EquilibriumDistribution"):
            self.fn("get", direction,
                    (u, f"Vector{D}<{T}>({T}{{0}})"),
                    (rho, f"{T}{{1}}"), returns=self.dtype)(
                *self.bind_vector(E.u_names, "u"),
                switch)

    # ---- namespace Equilibrium

    def gen_equilibrium(self):
        T = self.T
        E = self.equations
        u = self.vector("u")
        with self.sfg.namespace("Equilibrium"):
            self.fn("set", self.pdf_mut, u, self.rho_param, self.cell)(
                f"{T} &xyz0 = pdf_field->get(cell, {_ui(0)});",
                *self.bind_vector(E.u_names, "u"),
                *E.delta_rho.lines(),
                *(f"pdf_field->getF(&xyz0, {_ui(i)}) = "
                  f"{accessor_equations.print_expression(expr, self.dtype)};"
                  for i, expr in enumerate(E.equilibrium)))

    # ---- namespace Density

    def gen_density(self):
        T, D = self.T, self.D
        E = self.equations
        pdf, pdf_mut, cell, ci, values = (self.pdf_const, self.pdf_mut,
                                          self.cell, self.ci, self.values)
        density = self.density_param
        rho_in = self.sfg.var("rho_in", constify(self.dtype))
        scaled = E.density_scaled
        unshifted = E.unshifted_momentum_density
        used_u = E.used_pdfs(unshifted)
        used_s = E.used_pdfs(scaled)
        velocity_lines = [
            f"const {T} conversion = {T}{{1}} / rho;",
            f"Vector{D}<{T}> velocity;",
            *(f"velocity[{i}u] = md_{i} * conversion;" for i in range(D))]
        with self.sfg.namespace("Density"):
            self.fn("get", pdf, density, cell, returns=self.dtype)(
                f"const {T} &xyz0 = pdf_field->get(cell, {_ui(0)});",
                *self.pdf_loads(used_s),
                *scaled.lines(),
                "return rho;")

            self.fn("set", pdf_mut, rho_in, density, cell)(
                f"const {T} &xyz0 = pdf_field->get(cell, {_ui(0)});",
                *self.pdf_loads(used_u),
                *unshifted.lines(),
                "",
                "// calculate current velocity (before density change)",
                *velocity_lines,
                "",
                "Equilibrium::set(pdf_field, velocity, rho_in / density, cell);")

            self.fn("get", pdf, density, ci, returns=self.vector_return)(
                f"std::vector<{T}> out;",
                "out.reserve(ci.numCells());",
                *_interval_loops(
                    f"const {T} &xyz0 = pdf_field->get(x, y, z, {_ui(0)});",
                    *self.pdf_loads(used_s),
                    *scaled.lines(),
                    "out.emplace_back(rho);"),
                "return out;")

            self.fn("set", pdf_mut, values, density, ci)(
                "assert(uint_c(values.size()) == ci.numCells());",
                "auto values_it = values.begin();",
                f"auto const density_inv = {T}{{1}} / density;",
                *_interval_loops(
                    f"const {T} &xyz0 = pdf_field->get(x, y, z, {_ui(0)});",
                    *self.pdf_loads(used_u),
                    *unshifted.lines(),
                    "",
                    "// calculate current velocity (before density change)",
                    *velocity_lines,
                    "",
                    "Equilibrium::set(pdf_field, velocity, "
                    "*values_it * density_inv, Cell{x, y, z});",
                    "++values_it;"))

    # ---- namespace Velocity

    def gen_velocity(self):
        T, D = self.T, self.D
        E = self.equations
        pdf, pdf_mut, cell, ci, values = (self.pdf_const, self.pdf_mut,
                                          self.cell, self.ci, self.values)
        vel = self.vec_mut("velocity_field")
        frc = self.vec_const("force_field")
        u = self.vector("u")
        md = E.momentum_density
        used_md = E.used_pdfs(md)
        dens = E.density
        used_rho = E.used_pdfs(dens)
        setter = E.velocity_setter
        vector_args = ", ".join(f"md_{i} * rho_inv" for i in range(D))
        with self.sfg.namespace("Velocity"):
            self.fn("get", pdf, frc, cell, returns=AUTO)(
                f"const {T} &xyz0 = pdf_field->get(cell, {_ui(0)});",
                *self.pdf_loads(used_md),
                *self.cell_coordinates(),
                *self.force_loads(md),
                *md.lines(),
                f"const {T} rho_inv = {T}{{1}} / rho;",
                "",
                f"return Vector{D}<{T}>({vector_args});")

            self.fn("get", pdf, frc, ci, returns=AUTO)(
                f"std::vector<{T}> out;",
                f"out.reserve(ci.numCells() * uint_t({D}u));",
                *_interval_loops(
                    f"const {T} &xyz0 = pdf_field->get(x, y, z, {_ui(0)});",
                    *self.pdf_loads(used_md),
                    *self.force_loads(md),
                    *md.lines(),
                    f"const {T} rho_inv = {T}{{1}} / rho;",
                    *(f"out.emplace_back(md_{i} * rho_inv);"
                      for i in range(D))),
                "return out;")

            self.fn("set", pdf_mut, vel, frc, u, cell)(
                f"const {T} &xyz0 = pdf_field->get(cell, {_ui(0)});",
                *self.pdf_loads(used_rho),
                *dens.lines(),
                "",
                *self.cell_coordinates(),
                *self.force_loads(setter),
                *self.bind_vector(E.u_in_names, "u"),
                *setter.lines(),
                *(f"velocity_field->get(x, y, z, {_ui(i)}) = u[{i}u];"
                  for i in range(D)),
                "",
                self.equilibrium_call("cell"))

            self.fn("set", pdf_mut, vel, frc, values, ci)(
                f"assert(uint_c(values.size()) == ci.numCells() * "
                f"uint_t({D}u));",
                "auto u = values.data();",
                *_interval_loops(
                    f"{T} &pdf_xyz0 = pdf_field->get(x, y, z, {_ui(0)});",
                    f"{T} &vel_xyz0 = velocity_field->get(x, y, z, {_ui(0)});",
                    *self.pdf_loads(used_rho, ptr="pdf_xyz0"),
                    *dens.lines(),
                    "",
                    *self.force_loads(setter),
                    *self.bind_vector(E.u_in_names, "u"),
                    *setter.lines(),
                    *(f"velocity_field->getF(&vel_xyz0, {_ui(i)}) = u[{i}u];"
                      for i in range(D)),
                    "",
                    f"std::advance(u, {D});",
                    "",
                    self.equilibrium_call("Cell{x, y, z}")))

    # ---- namespace Force

    def gen_force(self):
        T, D = self.T, self.D
        E = self.equations
        pdf = self.pdf_const
        vel = self.vec_mut("velocity_field")
        frc = self.vec_mut("force_field")
        force = self.vector("force")
        cell, ci, values = self.cell, self.ci, self.values
        density = self.density_param
        eqs = E.momentum_density_force_setter
        used = E.used_pdfs(eqs)
        bind_force = self.bind_vector(E.force_names, "force")
        with self.sfg.namespace("Force"):
            self.fn("set", pdf, vel, frc, force, density, cell)(
                f"{T} const &pdf_xyz0 = pdf_field->get(cell, {_ui(0)});",
                f"{T} &vel_xyz0 = velocity_field->get(cell, {_ui(0)});",
                f"{T} &laf_xyz0 = force_field->get(cell, {_ui(0)});",
                *self.pdf_loads(used, ptr="pdf_xyz0"),
                "",
                *bind_force,
                *eqs.lines(),
                f"auto const rho_inv = {T}{{1}} / rho;",
                "",
                *(f"force_field->getF(&laf_xyz0, {_ui(i)}) = "
                  f"force[{i}u] / density;" for i in range(D)),
                "",
                *(f"velocity_field->getF(&vel_xyz0, {_ui(i)}) = "
                  f"md_{i} * rho_inv;" for i in range(D)))

            self.fn("set", pdf, vel, frc, values, density, ci)(
                f"assert(uint_c(values.size()) == ci.numCells() * "
                f"uint_t({D}u));",
                "auto force = values.data();",
                f"auto const density_inv = {T}{{1}} / density;",
                *_interval_loops(
                    f"{T} const &pdf_xyz0 = pdf_field->get(x, y, z, "
                    f"{_ui(0)});",
                    f"{T} &vel_xyz0 = velocity_field->get(x, y, z, "
                    f"{_ui(0)});",
                    f"{T} &laf_xyz0 = force_field->get(x, y, z, {_ui(0)});",
                    *self.pdf_loads(used, ptr="pdf_xyz0"),
                    "",
                    *bind_force,
                    *eqs.lines(),
                    f"auto const rho_inv = {T}{{1}} / rho;",
                    "",
                    *(f"force_field->getF(&laf_xyz0, {_ui(i)}) = "
                      f"force[{i}u] * density_inv;" for i in range(D)),
                    "",
                    *(f"velocity_field->getF(&vel_xyz0, {_ui(i)}) = "
                      f"md_{i} * rho_inv;" for i in range(D)),
                    "",
                    f"std::advance(force, {D});"))

    # ---- namespace MomentumDensity

    def gen_momentum_density(self):
        T, D = self.T, self.D
        E = self.equations
        frc = self.vec_const("force_field")
        # only the shifted momentum density is needed: this skips the
        # density and the populations it depends on
        eqs = E.momentum_density.pruned([f"md_{i}" for i in range(D)])
        with self.sfg.namespace("MomentumDensity"):
            self.fn("reduce", self.pdf_const, frc, self.density_param,
                    returns=AUTO)(
                f"Vector{D}<{T}> momentumDensity({T}{{0}});",
                *_interior_loops(
                    f"const {T} &xyz0 = pdf_field->get(x, y, z, {_ui(0)});",
                    *self.pdf_loads(E.used_pdfs(eqs)),
                    *self.force_loads(eqs),
                    *eqs.lines(),
                    "",
                    *(f"momentumDensity[{i}u] += md_{i} * density;"
                      for i in range(D))),
                "return momentumDensity;")

    # ---- namespace PressureTensor

    def gen_pressure_tensor(self):
        T, D = self.T, self.D
        E = self.equations
        pdf, cell, ci = self.pdf_const, self.cell, self.ci
        density = self.density_param
        eqs = E.second_moment
        used = E.used_pdfs(eqs)
        comps = range(D * D)
        with self.sfg.namespace("PressureTensor"):
            self.fn("get", pdf, density, cell, returns=AUTO)(
                "const int x = cell[0];",
                "const int y = cell[1];",
                "const int z = cell[2];",
                *self.pdf_offset_loads(used),
                *eqs.lines(),
                "",
                f"Matrix{D}<{T}> pressureTensor;",
                *(f"pressureTensor[{i}u] = p_{i};" for i in comps),
                "return pressureTensor * density;")

            self.fn("get", pdf, density, ci, returns=AUTO)(
                f"std::vector<{T}> out;",
                f"out.reserve(ci.numCells() * uint_t({D * D}u));",
                *_interval_loops(
                    *self.pdf_offset_loads(used),
                    *eqs.lines(),
                    "",
                    *(f"out.emplace_back(p_{i} * density);" for i in comps)),
                "return out;")

            self.fn("reduce", pdf, density, returns=AUTO)(
                f"Matrix{D}<{T}> pressureTensor({T}{{0}});",
                *_interior_loops(
                    *self.pdf_offset_loads(used),
                    *eqs.lines(),
                    "",
                    *(f"pressureTensor[{i}u] += p_{i};" for i in comps)),
                "return pressureTensor * density;")


def generate_field_accessors(sfg, ctx, lb_method):
    """
    Add the complete set of field accessors of ``lb_method`` to the header
    ``sfg``, which must be header-only (see ``KernelFileGenerator``).
    """
    _FieldAccessors(sfg, ctx, lb_method).generate()
