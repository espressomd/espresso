#
# Copyright (C) 2022 The ESPResSo project
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

import unittest as ut
import unittest_decorators as utx
import espressomd
import espressomd.EKSpecies
import numpy as np


@utx.skipIfMissingFeatures(["WALBERLA"])
class EKReaction(ut.TestCase):
    BOX_L = 11.
    AGRID = 1.0
    INITIAL_DENSITY = 1.0
    DIFFUSION_COEFFICIENT = 0.1
    TIME = 1000

    system = espressomd.System(box_l=[BOX_L, BOX_L, BOX_L])
    system.time_step = 1.0
    system.cell_system.skin = 0.4

    def tearDown(self) -> None:
        self.system.ekcontainer.clear()
        self.system.ekreactions.clear()

    def analytic_density_base(
            self, time: float, coeffs, rate_constant: float, init_density: float) -> float:
        """
        Calculates the base density of a species after a given time of a reaction.
        The reaction is defined via the stoichiometric coefficient of the edducts.
        To calculate the effective species density this base density needs to be multiplied by its stoichiometric coefficient.
        For the product density, one needs to substract this density from the init density.
        """
        order = sum(coeffs)
        factor = rate_constant
        for coeff in coeffs:
            factor *= coeff**coeff
        init_dens_factor = init_density ** (1 - order)
        return (init_dens_factor + (order - 1) *
                factor * time)**(1 / (1 - order))

    def test_reaction_single(self):
        self.detail_test_reaction(single_precision=True)

    def test_reaction_double(self):
        self.detail_test_reaction(single_precision=False)

    def detail_test_reaction(self, single_precision: bool):

        relative_precision: float = 1E-6 if single_precision else 1E-7

        lattice = espressomd.lb.LatticeWalberla(
            n_ghost_layers=1, agrid=self.AGRID)

        eksolver = espressomd.EKSpecies.EKNone(lattice=lattice)

        self.system.ekcontainer.tau = 1.0

        reaction_rate: float = 1e-5

        stoech_coeffs = [2.0, 1.0, 1.2, 2.2]
        product_coeff = 1.5
        educt_species = []
        reactants = []
        for coeff in stoech_coeffs:
            species = espressomd.EKSpecies.EKSpecies(
                lattice=lattice, density=coeff * self.INITIAL_DENSITY, kT=0.0,
                diffusion=self.DIFFUSION_COEFFICIENT, valency=0.0,
                advection=False, friction_coupling=False,
                ext_efield=[0., 0., 0.], single_precision=single_precision)
            self.system.ekcontainer.add(species)
            reactants.append(
                espressomd.EKSpecies.EKReactant(
                    ekspecies=species,
                    stoech_coeff=-coeff,
                    order=coeff))
            educt_species.append(species)

        ek_species_product = espressomd.EKSpecies.EKSpecies(
            lattice=lattice, density=0.0, diffusion=self.DIFFUSION_COEFFICIENT,
            kT=0.0, valency=0.0, advection=False, friction_coupling=False,
            ext_efield=[0., 0., 0.], single_precision=single_precision)
        self.system.ekcontainer.add(ek_species_product)
        reactants.append(
            espressomd.EKSpecies.EKReactant(
                ekspecies=ek_species_product,
                stoech_coeff=product_coeff,
                order=0.0))

        self.system.ekcontainer.solver = eksolver

        reaction = espressomd.EKSpecies.EKBulkReaction(
            reactants=reactants, coefficient=reaction_rate, lattice=lattice)

        self.system.ekreactions.add(reaction)

        self.system.integrator.run(self.TIME)

        domain_volume = np.product(ek_species_product.shape)

        measured_educt_densities = np.zeros(len(stoech_coeffs))
        for i, educt in enumerate(educt_species):
            measured_educt_densities[i] = np.sum(
                educt[:, :, :].density) / domain_volume
        measured_product_density = np.sum(
            ek_species_product[:, :, :].density) / domain_volume

        analytic_educt_densities = np.zeros(len(stoech_coeffs))
        for i, coeff in enumerate(stoech_coeffs):
            analytic_educt_densities[i] = coeff * self.analytic_density_base(
                self.TIME + 0.5, stoech_coeffs, reaction_rate, self.INITIAL_DENSITY)
        analytic_product_density = product_coeff * \
            (self.INITIAL_DENSITY -
                self.analytic_density_base(self.TIME + 0.5, stoech_coeffs,
                                           reaction_rate, self.INITIAL_DENSITY))

        np.testing.assert_allclose(
            measured_educt_densities / measured_educt_densities[0],
            analytic_educt_densities / analytic_educt_densities[0],
            rtol=relative_precision, atol=0)
        np.testing.assert_allclose(
            measured_educt_densities,
            analytic_educt_densities,
            rtol=2 * reaction_rate,
            atol=0)
        np.testing.assert_allclose(
            measured_product_density,
            analytic_product_density,
            rtol=2 *
            reaction_rate *
            len(stoech_coeffs),
            atol=0)


if __name__ == "__main__":
    ut.main()
