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
import sympy as sp
import numpy as np
import typing

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
    def __init__(self, dim, density_field, flux_field, diffusion, kT=None, velocity_field=None,
                 force_field=None, potential_field=None, valency=None, ext_efield=None):
        assert not ps.FieldType.is_staggered(density_field)

        if velocity_field is not None:
            assert not ps.FieldType.is_staggered(velocity_field)

        if force_field is not None:
            assert not ps.FieldType.is_staggered(force_field)

        if potential_field is not None:
            assert not ps.FieldType.is_staggered(potential_field)

        assert ps.FieldType.is_staggered(flux_field)

        self.dim = dim
        self.density_field = density_field
        self.velocity_field = velocity_field
        self.flux_field = flux_field
        self.diffusion = diffusion
        self.kT = kT
        self.force_field = force_field
        self.potential_field = potential_field
        self.valency = valency
        self.ext_efield = ext_efield

        full_stencil = ["C"] + self.flux_field.staggered_stencil + list(
            map(inverse_direction_string, self.flux_field.staggered_stencil))
        self.stencil = tuple(map(lambda d: tuple(
            ps.stencil.direction_string_to_offset(d, self.dim)), full_stencil))

        flux_expression = -self.diffusion * sp.Matrix(
            [ps.fd.diff(self.density_field, i) for i in range(self.density_field.spatial_dimensions)])

        if self.potential_field is not None and self.valency is not None:
            if ext_efield is not None:
                field = sp.Matrix([ps.fd.diff(self.potential_field, i) - ext_efield[i]
                                   for i in range(self.density_field.spatial_dimensions)])
            else:
                field = sp.Matrix([ps.fd.diff(self.potential_field, i)
                                   for i in range(self.density_field.spatial_dimensions)])

            flux_expression += - self.diffusion / self.kT * \
                self.density_field.center * self.valency * field

        self.disc = ps.fd.FVM1stOrder(
            self.density_field, flux=flux_expression, source=0)

        if self.velocity_field is not None:
            self.vof = ps.fd.VOF(
                self.flux_field,
                self.velocity_field,
                self.density_field)

    def flux_advection(self):
        if self.velocity_field is not None:
            return [ps.Assignment(j_adv.lhs, j_adv.lhs + j_adv.rhs)
                    for j_adv in self.vof]

    def flux(self, include_vof: bool = False,
             include_fluctuations: bool = False,
             rng_node: typing.Optional[ps.rng.RNGBase] = None):

        _flux_collection = ps.AssignmentCollection(
            [self.disc.discrete_flux(self.flux_field)])

        if include_fluctuations:
            if rng_node is None:
                raise ValueError(
                    "rng_node not provided but fluctuations requested")

            block_offsets = tuple(
                ps.TypedSymbol(
                    "block_offset_{}".format(i),
                    np.uint32) for i in range(
                    self.dim))

            rng_symbol_gen = random_symbol(_flux_collection.subexpressions,
                                           dim=self.dim,
                                           rng_node=rng_node,
                                           seed=ps.TypedSymbol(
                                               "seed", np.uint32),
                                           offsets=block_offsets)

            stencil = self.flux_field.staggered_stencil
            stencil_offsets = list(
                map(lambda d: ps.stencil.direction_string_to_offset(d), stencil))

            for i, (val, d, rng_symb) in enumerate(
                    zip(stencil, stencil_offsets, rng_symbol_gen)):
                assert _flux_collection.main_assignments[i].lhs == self.flux_field.staggered_access(
                    val)
                _flux_collection.main_assignments[i] = ps.Assignment(
                    self.flux_field.staggered_access(val),
                    _flux_collection.main_assignments[i].rhs + sp.sqrt(
                        2 * self.diffusion * discretize(self.density_field.center, d)) / sp.Matrix(
                        d).norm() * rng_symb * sp.sqrt(
                        3) / 4)

        if include_vof:
            assert self.velocity_field is not None, "velocity field is not provided!"

            for i, j_adv in enumerate(self.vof):
                assert _flux_collection.main_assignments[i].lhs == j_adv.lhs
                _flux_collection.main_assignments[i] = ps.Assignment(
                    j_adv.lhs,
                    _flux_collection.main_assignments[i].rhs + j_adv.rhs)

        return _flux_collection

    def continuity(self):
        return self.disc.discrete_continuity(self.flux_field)

    def friction_coupling(self):
        if self.kT is None or self.force_field is None:
            raise RuntimeError("kT or f is not provided!")

        stencil = self.flux_field.staggered_stencil + \
            [ps.stencil.inverse_direction_string(
                d) for d in self.flux_field.staggered_stencil]

        return ps.AssignmentCollection([ps.Assignment(self.force_field.center_vector, self.kT / (2 * self.diffusion) * sum([self.flux_field.staggered_access(val) * sp.Matrix(
            ps.stencil.direction_string_to_offset(val)) for val in stencil[1:]], self.flux_field.staggered_access(stencil[0]) * sp.Matrix(ps.stencil.direction_string_to_offset(stencil[0]))))])


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
