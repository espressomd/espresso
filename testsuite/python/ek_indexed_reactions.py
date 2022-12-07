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
    BOX_L = [22., 2., 2.]
    PADDING = 1
    WIDTH = 20.
    AGRID = 1.0
    INITIAL_DENSITIES = [1.7, 1.3]
    DIFFUSION_COEFFICIENTS = np.array([0.4, 0.2])
    REACTION_RATES = np.array([4e-3, 8e-3])
    TIME = 6000

    system = espressomd.System(box_l=BOX_L)
    system.time_step = 1.0
    system.cell_system.skin = 0.4

    def tearDown(self) -> None:
        self.system.ekcontainer.clear()
        self.system.ekreactions.clear()

    def analytic_density_profiles(
            self, width, reaction_rates, diffusion_coefficients, initial_densities, agrid):
        rezipr_diff = 1 / \
            diffusion_coefficients[0] + 1 / diffusion_coefficients[1]
        rezipr_rate = 1 / reaction_rates[0] + 1 / reaction_rates[1]
        actual_width = width - agrid
        slopes = sum(initial_densities) / (diffusion_coefficients *
                                           (rezipr_rate + actual_width / 2 * rezipr_diff))
        midvalues = sum(initial_densities) / (reaction_rates * (rezipr_rate +
                                                                actual_width / 2 * rezipr_diff)) + actual_width / 2 * slopes

        x = np.linspace(-actual_width / 2, actual_width /
                        2, int(width / agrid))
        values_a = slopes[0] * x + midvalues[0]
        values_b = -slopes[1] * x + midvalues[1]
        return values_a, values_b

    def test_reaction_single(self):
        self.detail_test_reaction(single_precision=True)

    def test_reaction_double(self):
        self.detail_test_reaction(single_precision=False)

    def detail_test_reaction(self, single_precision: bool):

        lattice = espressomd.lb.LatticeWalberla(
            n_ghost_layers=1, agrid=self.AGRID)

        eksolver = espressomd.EKSpecies.EKNone(lattice=lattice)

        self.system.ekcontainer.tau = 1.0

        self.system.ekcontainer.solver = eksolver

        species_A = espressomd.EKSpecies.EKSpecies(
            lattice=lattice, density=self.INITIAL_DENSITIES[0], kT=0.0,
            diffusion=self.DIFFUSION_COEFFICIENTS[0], valency=0.0,
            advection=False, friction_coupling=False, ext_efield=[0, 0, 0],
            single_precision=single_precision)
        self.system.ekcontainer.add(species_A)

        species_B = espressomd.EKSpecies.EKSpecies(
            lattice=lattice, density=self.INITIAL_DENSITIES[1], kT=0.0,
            diffusion=self.DIFFUSION_COEFFICIENTS[1], valency=0.0,
            advection=False, friction_coupling=False, ext_efield=[0, 0, 0],
            single_precision=single_precision)
        self.system.ekcontainer.add(species_B)

        coeffs_left = [-1.0, 1.0]
        reactants_left = []
        reactants_left.append(
            espressomd.EKSpecies.EKReactant(
                ekspecies=species_A,
                stoech_coeff=coeffs_left[0],
                order=1.0))
        reactants_left.append(
            espressomd.EKSpecies.EKReactant(
                ekspecies=species_B,
                stoech_coeff=coeffs_left[1],
                order=0.0))

        reaction_left = espressomd.EKSpecies.EKIndexedReaction(
            reactants=reactants_left, coefficient=self.REACTION_RATES[0],
            lattice=lattice)
        reaction_left[1, :, :] = True

        coeffs_right = [1.0, -1.0]
        reactants_right = []
        reactants_right.append(
            espressomd.EKSpecies.EKReactant(
                ekspecies=species_A,
                stoech_coeff=coeffs_right[0],
                order=0.0))
        reactants_right.append(
            espressomd.EKSpecies.EKReactant(
                ekspecies=species_B,
                stoech_coeff=coeffs_right[1],
                order=1.0))

        reaction_right = espressomd.EKSpecies.EKIndexedReaction(
            reactants=reactants_right, coefficient=self.REACTION_RATES[1],
            lattice=lattice)
        reaction_right[-2, :, :] = True

        self.system.ekreactions.add(reaction_left)
        self.system.ekreactions.add(reaction_right)

        wall_left = espressomd.shapes.Wall(normal=[1, 0, 0], dist=self.PADDING)
        wall_right = espressomd.shapes.Wall(
            normal=[-1, 0, 0], dist=-self.BOX_L[0] + self.PADDING)
        for obj in (wall_left, wall_right):
            species_A.add_boundary_from_shape(
                obj, [0.0, 0.0, 0.0], espressomd.EKSpecies.FluxBoundary)
            species_B.add_boundary_from_shape(
                obj, [0.0, 0.0, 0.0], espressomd.EKSpecies.FluxBoundary)

        self.system.integrator.run(self.TIME)

        density_profile = np.zeros((2, int(self.WIDTH / self.AGRID)))
        for x in range(int(self.WIDTH / self.AGRID)):
            density_profile[0, x] = np.mean(
                species_A[x + self.PADDING, :, :].density)
            density_profile[1, x] = np.mean(
                species_B[x + self.PADDING, :, :].density)

        analytic_density_profile = np.zeros((2, int(self.WIDTH / self.AGRID)))
        analytic_density_profile[0], analytic_density_profile[1] = \
            self.analytic_density_profiles(self.WIDTH,
                                           self.REACTION_RATES,
                                           self.DIFFUSION_COEFFICIENTS,
                                           self.INITIAL_DENSITIES,
                                           self.AGRID)

        np.testing.assert_allclose(
            density_profile,
            analytic_density_profile,
            rtol=self.REACTION_RATES[0],
            atol=0)


if __name__ == "__main__":
    ut.main()
