#
# Copyright (C) 2026 The ESPResSo project
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
Equations of the LB field accessors, derived from an lbmpy method as sympy
expressions; `print_expression` prints an expression as C++ with
pystencils 2.1.

All equations are written with plain symbols; the code that embeds them
must replace these symbols (see `lb_accessors.LbAccessors.reads`):

* ``f_0 ... f_{Q-1}``: populations of the cell,
* ``force_0 ... force_{D-1}``: components of the force density
  (replaces the force-model field accesses ``force[0,0,0](i)``),
* ``u_in_0 ... u_in_{D-1}``: velocity passed by the caller,
* ``rho`` / ``delta_rho`` / ``u_0 ... u_{D-1}``: density and velocity
  inputs of the equilibrium,
* ``density``: conversion factor between lattice and MD units.

The equations introduce temporary names such as ``vel0Term``,
``momdensity_i``, ``md_i`` and ``p_i``; the embedding code must not reuse
them for other variables.
"""

import dataclasses
import functools

import pystencils as ps
import sympy as sp
# NOTE: the ``pystencils.backend.*`` modules are pystencils-internal APIs,
# tied to the tested pystencils version (see ``TESTED_TOOLCHAIN`` in
# ``sfg_driver.py``). Re-check them whenever the toolchain is upgraded.
from pystencils.backend.emission import CAstPrinter
from pystencils.backend.kernelcreation import AstFactory, KernelCreationContext
from pystencils.backend.transformations import EliminateConstants
from pystencils.stencil import offset_to_direction_string

__all__ = ["EquationSet", "AccessorEquations", "print_expression"]


def _symbol_names(expr):
    return {s.name for s in expr.free_symbols}


def _printing_context(dtype):
    ctx = KernelCreationContext(default_dtype=dtype)
    return (AstFactory(ctx), EliminateConstants(ctx, fold_floats=True),
            CAstPrinter())


def print_expression(expr, dtype):
    """Print a sympy expression as a C++ expression typed with ``dtype``.

    The printer types every plain symbol with the working float type; there
    are no integer symbols.
    """
    factory, fold, printer = _printing_context(dtype)
    return printer(fold(factory.parse_sympy(sp.sympify(expr))))


@dataclasses.dataclass(frozen=True)
class EquationSet:
    """
    Ordered list of assignments to plain symbols.

    Attributes:
        assignments: tuple of ``(lhs symbol, rhs expression)``
    """
    assignments: tuple

    @property
    def outputs(self):
        return [lhs.name for lhs, _ in self.assignments]

    @property
    def inputs(self):
        """Names of the symbols that are read but never assigned."""
        inputs, defined = set(), set()
        for lhs, rhs in self.assignments:
            inputs |= _symbol_names(rhs) - defined
            defined.add(lhs.name)
        return inputs

    def pruned(self, targets):
        """
        Keep only the assignments needed to compute the symbols ``targets``
        (names). The relative order of the assignments is preserved.
        """
        targets = set(targets)
        unknown = targets - set(self.outputs)
        if unknown:
            raise ValueError(
                f"cannot prune to unknown targets {sorted(unknown)}; "
                f"available outputs: {self.outputs}")
        needed = set(targets)
        kept = []
        for lhs, rhs in reversed(self.assignments):
            if lhs.name in needed:
                kept.append((lhs, rhs))
                needed.discard(lhs.name)
                needed |= _symbol_names(rhs)
        return dataclasses.replace(self, assignments=tuple(reversed(kept)))


class AccessorEquations:
    """
    All equations of the field accessors of an lbmpy method.

    Args:
        lb_method: lbmpy method with a force model whose force vector is a
            tuple of field accesses; zero-centered PDFs and a compressible
            fluid are required
    """

    def __init__(self, lb_method):
        if lb_method.force_model is None:
            raise NotImplementedError(
                "the accessors require a force model")
        cqc = lb_method.conserved_quantity_computation
        if not (cqc.compressible and cqc.zero_centered_pdfs):
            raise NotImplementedError(
                "the accessors support compressible fluids with "
                "zero-centered PDFs only")
        self.lb_method = lb_method
        self.D = lb_method.stencil.D
        self.Q = lb_method.stencil.Q
        self.stencil = lb_method.stencil
        self._cqc = cqc

        self.pdf_names = [f"f_{i}" for i in range(self.Q)]
        self.force_names = [f"force_{i}" for i in range(self.D)]
        self.u_in_names = [f"u_in_{i}" for i in range(self.D)]
        self.u_names = [f"u_{i}" for i in range(self.D)]

        self._pdfs = sp.symbols(self.pdf_names)
        self._rho = sp.Symbol("rho")
        self._delta_rho = sp.Symbol("delta_rho")
        self._density = sp.Symbol("density")
        self._md = sp.symbols([f"md_{i}" for i in range(self.D)])
        self._p = sp.symbols([f"p_{i}" for i in range(self.D**2)])
        self._u = sp.symbols(self.u_names)
        self._u_in = sp.symbols(self.u_in_names)

        # replace the force-model field accesses by plain symbols
        force_vector = list(lb_method.force_model.symbolic_force_vector)
        if len(force_vector) != self.D:
            raise ValueError(
                f"the force vector has {len(force_vector)} components, "
                f"expected {self.D}")
        if not all(isinstance(f, ps.Field.Access) for f in force_vector):
            raise NotImplementedError(
                "the force vector must consist of field accesses")
        self._force = sp.symbols(self.force_names)
        self._force_subs = dict(zip(force_vector, self._force))

    # -- helpers

    def _equation_set(self, assignments, subs=None):
        pairs = []
        for assignment in assignments:
            lhs, rhs = assignment.lhs, assignment.rhs
            rhs = rhs.xreplace(self._force_subs)
            if subs:
                rhs = rhs.xreplace(subs)
            pairs.append((lhs, rhs))
        return EquationSet(tuple(pairs))

    def _output_equations(self, outputs):
        ac = self._cqc.output_equations_from_pdfs(self._pdfs, outputs)
        return self._equation_set(ac.all_assignments)

    @staticmethod
    def _replace_rhs(eqs, name, function):
        return dataclasses.replace(eqs, assignments=tuple(
            (lhs, function(rhs)) if lhs.name == name else (lhs, rhs)
            for lhs, rhs in eqs.assignments))

    # -- equation sets used by the accessors

    @functools.cached_property
    def density(self):
        """``rho`` in lattice units."""
        return self._output_equations({"density": self._rho})

    @functools.cached_property
    def density_scaled(self):
        """``rho`` multiplied by the input symbol ``density``."""
        return self._replace_rhs(self.density, "rho",
                                 lambda rhs: rhs * self._density)

    @functools.cached_property
    def momentum_density(self):
        """Half-force-shifted momentum density ``md_i`` and ``rho``."""
        return self._output_equations(
            {"density": self._rho, "momentum_density": self._md})

    @functools.cached_property
    def momentum_density_force_setter(self):
        """Same as `momentum_density`, with ``force_i`` divided by
        ``density``, since the force is given in MD units."""
        subs = {f: f / self._density for f in self._force}
        eqs = self.momentum_density
        return dataclasses.replace(eqs, assignments=tuple(
            (lhs, rhs.xreplace(subs)) for lhs, rhs in eqs.assignments))

    @functools.cached_property
    def unshifted_momentum_density(self):
        """Momentum density without force shift (``md_i``) and ``rho``."""
        eqs = self.momentum_density
        zero = {f: 0 for f in self._force}
        return dataclasses.replace(eqs, assignments=tuple(
            (lhs, rhs.xreplace(zero)) for lhs, rhs in eqs.assignments)
        ).pruned(["rho"] + [m.name for m in self._md])

    @functools.cached_property
    def second_moment(self):
        """Second moment ``p_i`` (``D x D``, row major) of the PDFs."""
        return self._output_equations({"moment2": self._p})

    @functools.cached_property
    def velocity_setter(self):
        """
        Velocity ``u_i`` of the equilibrium, from the velocity ``u_in_i``
        requested by the caller and the density ``rho``: the requested
        velocity is shifted by half the force.
        """
        eqin = self._cqc.equilibrium_input_equations_from_init_values(
            self._rho, self._u)
        by_name = {a.lhs.name: a for a in eqin.all_assignments}
        # the lbmpy velocity symbols are both the inputs and the outputs
        subs = dict(zip(self._u, self._u_in))
        pairs = []
        for u in self._u:
            assignment = by_name[u.name]
            pairs.append((u, assignment.rhs.xreplace(self._force_subs)
                          .xreplace(subs)))
        return EquationSet(tuple(pairs))

    def _equilibrium_expressions(self):
        eq = self.lb_method.get_equilibrium()
        exprs = [a.rhs for a in eq.main_assignments]
        allowed = {"rho", "delta_rho", *self.u_names}
        for expr in exprs:
            if not _symbol_names(expr) <= allowed:
                raise NotImplementedError(
                    "the equilibrium depends on unsupported symbols: "
                    f"{sorted(_symbol_names(expr) - allowed)}")
        return exprs

    @functools.cached_property
    def equilibrium(self):
        """Expressions of the equilibrium populations, in terms of
        ``rho``, ``delta_rho`` and ``u_i``."""
        return self._equilibrium_expressions()

    @functools.cached_property
    def equilibrium_from_direction(self):
        """
        Pairs of the waLBerla direction name and the full (not delta)
        equilibrium, in terms of ``rho`` and ``u_i``: ``delta_rho`` is
        replaced by ``rho``. This assumes that lbmpy's ``delta_rho`` is
        ``rho - 1`` (zero-centered PDFs); ``NumericTest`` guards it.
        """
        exprs = [e.xreplace({self._delta_rho: self._rho})
                 for e in self._equilibrium_expressions()]
        names = [offset_to_direction_string(d) for d in self.stencil]
        return list(zip(names, exprs))

    # -- stencil

    @functools.cached_property
    def inverse_neighbor_offsets(self):
        """Offsets ``-c_i`` of the neighbor cell the population ``i`` is
        pulled from."""
        return [tuple(-c for c in direction) for direction in self.stencil]

    def used_pdfs(self, eqs):
        """Indices of the populations an equation set reads."""
        names = eqs.inputs
        return [i for i, name in enumerate(self.pdf_names) if name in names]

    def used_forces(self, eqs):
        """Indices of the force components an equation set reads."""
        names = eqs.inputs
        return [i for i, name in enumerate(self.force_names) if name in names]
