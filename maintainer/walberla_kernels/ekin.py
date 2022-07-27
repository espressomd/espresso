import pystencils as ps
import sympy as sp
import numpy as np
from pystencils.fd.derivation import FiniteDifferenceStaggeredStencilDerivation
from pystencils.fd.finitevolumes import get_access_and_direction

from pystencils.rng import random_symbol

from pystencils.stencil import inverse_direction_string


# this is from ps.fd.FVM1stOrder.discrete_flux.discretize
def discretize(term, neighbor):
    if isinstance(term, sp.Matrix):
        nw = term.applyfunc(lambda t: discretize(t, neighbor))
        return nw
    elif isinstance(term, ps.field.Field.Access):
        avg = (term.get_shifted(*neighbor) + term) * sp.Rational(1, 2)
        return avg
    elif isinstance(term, ps.fd.Diff):
        access, direction = get_access_and_direction(term)

        fds = FiniteDifferenceStaggeredStencilDerivation(
            neighbor, access.field.spatial_dimensions, direction)
        return fds.apply(access)

    if term.args:
        new_args = [discretize(a, neighbor) for a in term.args]
        return term.func(*new_args)
    else:
        return term


class EK:
    def __init__(self, dim, rho, j, diffusion, kT=None, u=None,
                 f=None, phi=None, valency=None, ext_efield=None):
        assert not ps.FieldType.is_staggered(rho)

        if u is not None:
            assert not ps.FieldType.is_staggered(u)

        if f is not None:
            assert not ps.FieldType.is_staggered(f)

        if phi is not None:
            assert not ps.FieldType.is_staggered(phi)

        assert ps.FieldType.is_staggered(j)

        self.dim = dim
        self.rho = rho
        self.u = u
        self.j = j
        self.diffusion = diffusion
        self.kT = kT
        self.f = f
        self.phi = phi
        self.valency = valency
        self.ext_efield = ext_efield

        full_stencil = ["C"] + self.j.staggered_stencil + \
            list(map(inverse_direction_string, self.j.staggered_stencil))
        self.stencil = tuple(map(lambda d: tuple(
            ps.stencil.direction_string_to_offset(d, self.dim)), full_stencil))

        flux_expression = -self.diffusion * sp.Matrix(
            [ps.fd.diff(self.rho, i) for i in range(self.rho.spatial_dimensions)])

        if self.phi is not None and self.valency is not None:
            if ext_efield is not None:
                field = sp.Matrix([ps.fd.diff(self.phi, i) - ext_efield[i]
                                   for i in range(self.rho.spatial_dimensions)])
            else:
                field = sp.Matrix([ps.fd.diff(self.phi, i)
                                   for i in range(self.rho.spatial_dimensions)])

            flux_expression += - self.diffusion / self.kT * \
                self.rho.center * self.valency * field

        self.disc = ps.fd.FVM1stOrder(self.rho, flux=flux_expression, source=0)

        if self.u is not None:
            self.vof = ps.fd.VOF(self.j, self.u, self.rho)

    def flux_advection(self):
        if self.u is not None:
            return [ps.Assignment(j_adv.lhs, j_adv.lhs + j_adv.rhs)
                    for j_adv in self.vof]

    def flux(self, include_vof: bool = False,
             include_fluctuations: bool = False):

        _flux_collection = ps.AssignmentCollection(
            [self.disc.discrete_flux(self.j)])

        if include_fluctuations:

            rng_symbol_gen = random_symbol(_flux_collection.subexpressions,
                                           dim=self.dim,
                                           rng_node=ps.rng.PhiloxTwoDoubles,
                                           seed=ps.TypedSymbol("seed", np.uint32))

            stencil = self.j.staggered_stencil
            stencil_offsets = list(
                map(lambda d: ps.stencil.direction_string_to_offset(d), stencil))

            for i, (val, d, rng_symb) in enumerate(
                    zip(stencil, stencil_offsets, rng_symbol_gen)):
                assert _flux_collection.main_assignments[i].lhs == self.j.staggered_access(
                    val)
                _flux_collection.main_assignments[i] = ps.Assignment(
                    self.j.staggered_access(val),
                    _flux_collection.main_assignments[i].rhs + sp.sqrt(
                        2 * self.diffusion * discretize(self.rho.center, d)) / sp.Matrix(d).norm() * rng_symb * sp.sqrt(
                        3) / 4)

        if include_vof:
            assert self.u is not None, "velocity field is not provided!"

            for i, j_adv in enumerate(self.vof):
                assert _flux_collection.main_assignments[i].lhs == j_adv.lhs
                _flux_collection.main_assignments[i] = ps.Assignment(
                    j_adv.lhs,
                    _flux_collection.main_assignments[i].rhs + j_adv.rhs)

        return _flux_collection

    def continuity(self):
        return self.disc.discrete_continuity(self.j)

    def friction_coupling(self):
        if self.kT is None or self.f is None:
            raise RuntimeError("kT or f is not provided!")

        stencil = self.j.staggered_stencil + \
            [ps.stencil.inverse_direction_string(
                d) for d in self.j.staggered_stencil]

        return ps.Assignment(self.f.center_vector, self.kT / (2 * self.diffusion) * sum(
            [self.j.staggered_access(val) * sp.Matrix(ps.stencil.direction_string_to_offset(val)) for val in
             stencil[1:]],
            self.j.staggered_access(stencil[0]) * sp.Matrix(ps.stencil.direction_string_to_offset(stencil[0]))))


class Reaction:
    def __init__(self, species, orders, stoechom_coefs, rate_coef):
        self.species = species
        self.orders = orders
        self.stoechom_coefs = stoechom_coefs
        self.rate_coef = rate_coef

    def generate_reaction(self, num_reactants: int) -> ps.AssignmentCollection:
        if num_reactants > len(self.species):
            raise ValueError(
                "Not enough species defined for number of requested reactants")

        # read density fields into subexpressions
        rho_symbols = sp.symbols(f"local_rho_:{num_reactants}")
        rate_symbol = sp.Symbol("rate_factor")

        subexpressions = [
            ps.Assignment(
                rho_symbols[i],
                self.species[i].center) for i in range(num_reactants)]

        rate = self.rate_coef
        for i in range(num_reactants):
            rate *= sp.Pow(rho_symbols[i], self.orders[i])

        subexpressions.append(ps.Assignment(rate_symbol, rate))

        main_assignments = []
        for i in range(num_reactants):
            main_assignments.append(ps.Assignment(self.species[i].center,
                                                  rho_symbols[i] + rate_symbol * self.stoechom_coefs[i]))

        collection = ps.AssignmentCollection(subexpressions=subexpressions,
                                             main_assignments=main_assignments)

        return collection
