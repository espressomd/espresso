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

import unittest as ut
import unittest_decorators as utx
import espressomd
import espressomd.electrokinetics
import numpy as np


class EKTest:
    BOX_L = 11.
    AGRID = 1.1
    INITIAL_DENSITY = 1.0
    DIFFUSION_COEFFICIENT = 0.1
    TIME = 500
    TAU = 1.9

    system = espressomd.System(box_l=[BOX_L, BOX_L, BOX_L])
    system.time_step = TAU
    system.cell_system.skin = 0.4

    def tearDown(self):
        self.system.ekcontainer = None

    def analytic_density_base(
            self, time: float, coeffs, rate_constant: float, init_density: float) -> float:
        """
        Calculates the base density of a species after a given time of a reaction.
        The reaction is defined via the stoichiometric coefficient of the educts.
        To calculate the effective species density this base density needs to be multiplied by its stoichiometric coefficient.
        For the product density, one needs to subtract this density from the init density.
        """
        order = sum(coeffs)
        factor = rate_constant
        for coeff in coeffs:
            factor *= coeff**coeff
        init_dens_factor = init_density ** (1 - order)
        return (init_dens_factor + (order - 1) *
                factor * time)**(1 / (1 - order))

    def test_reaction(self):

        relative_precision = 1E-6 if self.ek_params["single_precision"] else 1E-7

        lattice = espressomd.electrokinetics.Lattice(
            n_ghost_layers=1, agrid=self.AGRID)

        eksolver = espressomd.electrokinetics.EKNone(
            lattice=lattice, tau=self.TAU)
        self.system.ekcontainer = espressomd.electrokinetics.EKContainer(
            tau=self.TAU, solver=eksolver)

        reaction_rate: float = 1e-5

        stoech_coeffs = [2.0, 1.0, 1.2, 2.2]
        product_coeff = 1.5
        educt_species = []
        reactants = []
        for coeff in stoech_coeffs:
            species = self.ek_species_class(
                lattice=lattice, density=coeff * self.INITIAL_DENSITY,
                diffusion=self.DIFFUSION_COEFFICIENT, valency=0.0,
                advection=False, friction_coupling=False, tau=self.TAU,
                **self.ek_params)
            self.system.ekcontainer.add(species)
            reactants.append(
                espressomd.electrokinetics.EKReactant(
                    ekspecies=species,
                    stoech_coeff=-coeff,
                    order=coeff))
            educt_species.append(species)

        ek_species_product = self.ek_species_class(
            lattice=lattice, density=0.0, diffusion=self.DIFFUSION_COEFFICIENT,
            valency=0.0, advection=False, friction_coupling=False,
            tau=self.TAU, **self.ek_params)
        self.system.ekcontainer.add(ek_species_product)
        reactants.append(
            espressomd.electrokinetics.EKReactant(
                ekspecies=ek_species_product,
                stoech_coeff=product_coeff,
                order=0.0))

        reaction = espressomd.electrokinetics.EKBulkReaction(
            reactants=reactants, coefficient=reaction_rate, lattice=lattice, tau=self.TAU)

        self.system.ekcontainer.reactions.add(reaction)
        self.assertEqual(len(self.system.ekcontainer.reactions), 1)

        self.system.integrator.run(self.TIME)

        domain_volume = np.prod(np.copy(ek_species_product.shape))
        analytic_time = (self.TIME + 0.5) * self.system.time_step

        measured_educt_densities = np.zeros(len(stoech_coeffs))
        for i, educt in enumerate(educt_species):
            measured_educt_densities[i] = np.copy(np.sum(
                educt[:, :, :].density)) / domain_volume
        measured_product_density = np.copy(np.sum(
            ek_species_product[:, :, :].density)) / domain_volume

        analytic_educt_densities = np.zeros(len(stoech_coeffs))
        for i, coeff in enumerate(stoech_coeffs):
            analytic_educt_densities[i] = coeff * self.analytic_density_base(
                analytic_time, stoech_coeffs, reaction_rate, self.INITIAL_DENSITY)
        analytic_product_density = product_coeff * \
            (self.INITIAL_DENSITY -
                self.analytic_density_base(analytic_time, stoech_coeffs,
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

        self.system.ekcontainer.reactions.remove(reaction)
        self.assertEqual(len(self.system.ekcontainer.reactions), 0)


@utx.skipIfMissingFeatures(["WALBERLA"])
class EKReactionDoublePrecisionCPU(EKTest, ut.TestCase):
    ek_species_class = espressomd.electrokinetics.EKSpecies
    ek_params = {"single_precision": False, "gpu": False}


@utx.skipIfMissingFeatures(["WALBERLA"])
class EKReactionSinglePrecisionCPU(EKTest, ut.TestCase):
    ek_species_class = espressomd.electrokinetics.EKSpecies
    ek_params = {"single_precision": True, "gpu": False}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class EKReactionDoublePrecisionGPU(EKTest, ut.TestCase):
    ek_species_class = espressomd.electrokinetics.EKSpecies
    ek_params = {"single_precision": False, "gpu": True}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class EKReactionSinglePrecisionGPU(EKTest, ut.TestCase):
    ek_species_class = espressomd.electrokinetics.EKSpecies
    ek_params = {"single_precision": True, "gpu": True}


if __name__ == "__main__":
    ut.main()
