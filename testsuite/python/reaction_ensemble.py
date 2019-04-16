#
# Copyright (C) 2013-2018 The ESPResSo project
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

"""Testmodule for the Reaction Ensemble.
"""
import os
import sys
import unittest as ut
import numpy as np
import espressomd  # pylint: disable=import-error
from espressomd import reaction_ensemble


class ReactionEnsembleTest(ut.TestCase):

    """Test the core implementation of the reaction ensemble."""

    # The reaction ensemble follows the ideal titration curve only if N>>1,
    # Ideal curve is derived in the grandcanionical ensemble and for low N
    # there are systematic devations caused by differences between the
    # ensembles. This is not an error but a fundamental difference (various
    # ensembles are equivalent only in the thermodynamic limit N \to \infty
    N0 = 40
    c0 = 0.00028
    type_HA = 0
    type_A = 1
    type_H = 2
    target_alpha = 0.6
    # We get best statistics at alpha=0.5 Then the test is least sensistive to # the exact sequence of random numbers and does not require hard-coded
    # output values
    temperature = 1.0
    exclusion_radius = 1.0
    # could be in this test for example anywhere in the range 0.000001 ... 9,
    reactant_types = [type_HA]
    reactant_coefficients = [1]
    product_types = [type_A, type_H]
    product_coefficients = [1, 1]
    nubar = 1
    system = espressomd.System(box_l=np.ones(3) * (N0 / c0)**(1.0 / 3.0))
    system.seed = system.cell_system.get_state()['n_nodes'] * [67]
    np.random.seed(69)  # make reaction code fully deterministic
    system.cell_system.skin = 0.4
    volume = np.prod(system.box_l)  # cuboid box
    # Calculate gamma which should lead to target_alpha with given N0 and V
    # Requires N0>>1, otherwise discrete character of N changes the statistics (N>20 should suffice)
    # gamma = prod_i (N_i / V) = alpha^2 N0 / (1-alpha)*V**(-nubar)
    # degree of dissociation alpha = N_A / N_HA = N_H / N_0
    gamma = target_alpha**2 / (1. - target_alpha) * N0 / (volume**nubar)
    RE = reaction_ensemble.ReactionEnsemble(
        temperature=temperature,
        exclusion_radius=exclusion_radius, seed=12)

    @classmethod
    def setUpClass(cls):
        """Prepare a testsystem."""
        for i in range(0, 2 * cls.N0, 2):
            cls.system.part.add(id=i, pos=np.random.random(
                3) * cls.system.box_l, type=cls.type_A)
            cls.system.part.add(id=i + 1, pos=np.random.random(3) *
                                cls.system.box_l, type=cls.type_H)

        cls.RE.add_reaction(
            gamma=cls.gamma,
            reactant_types=cls.reactant_types,
            reactant_coefficients=cls.reactant_coefficients,
            product_types=cls.product_types,
            product_coefficients=cls.product_coefficients,
            default_charges={cls.type_HA: 0, cls.type_A: -1, cls.type_H: +1}, check_for_electroneutrality=True)

    @classmethod
    def ideal_alpha(cls, gamma, N0, V, nubar):
        # gamma = prod_i (N_i / V) = alpha^2 N0 / (1-alpha)*V**(-nubar)
        # degree of dissociation alpha = N_A / N_HA = N_H / N_0
        X = 2 * N0 / (gamma * V**nubar)
        return (np.sqrt(1 + 2 * X) - 1) / X

    def test_ideal_titration_curve(self):
        N0 = ReactionEnsembleTest.N0
        type_A = ReactionEnsembleTest.type_A
        type_H = ReactionEnsembleTest.type_H
        type_HA = ReactionEnsembleTest.type_HA
        box_l = ReactionEnsembleTest.system.box_l
        system = ReactionEnsembleTest.system
        gamma = ReactionEnsembleTest.gamma
        nubar = ReactionEnsembleTest.nubar

        volume = ReactionEnsembleTest.volume
        RE = ReactionEnsembleTest.RE
        target_alpha = ReactionEnsembleTest.target_alpha

        # chemical warmup - get close to chemical equilibrium before we start
        # sampling
        RE.reaction(20 * N0)

        average_NH = 0.0
        average_NHA = 0.0
        average_NA = 0.0
        num_samples = 1000
        for i in range(num_samples):
            RE.reaction(10)
            average_NH += system.number_of_particles(type=type_H)
            average_NHA += system.number_of_particles(type=type_HA)
            average_NA += system.number_of_particles(type=type_A)
        average_NH /= num_samples
        average_NA /= num_samples
        average_NHA /= num_samples
        average_alpha = average_NA / float(N0)
        print(average_alpha)
        # Note: with 40 particles, alpha=0.5 and 1000*10 reactions, standard
        # deviation of average alpha is about 0.003 (determined from 40
        # repeated simulations).  We set the desired accuracy to 5*std = 0.015
        rel_error_alpha = abs(
            average_alpha - target_alpha) / target_alpha
        # relative error
        self.assertLess(
            rel_error_alpha,
            0.015,
            msg="\nDeviation from ideal titration curve is too big for the given input parameters.\n"
            + "  gamma: " + str(gamma)
            + "  average_NH: " + str(average_NH)
            + "  average_NA: " + str(average_NA)
            + "  average_NHA:" + str(average_NHA)
            + "  average alpha: " + str(average_alpha)
            + "  target_alpha: " + str(target_alpha)
        )

    def test_reaction_system(self):
        RE_status = ReactionEnsembleTest.RE.get_status()
        forward_reaction = RE_status["reactions"][0]
        for i in range(len(forward_reaction["reactant_types"])):
            self.assertEqual(
                ReactionEnsembleTest.reactant_types[i],
                forward_reaction["reactant_types"][i],
                msg="reactant type not set correctly.")
        for i in range(len(forward_reaction["reactant_coefficients"])):
            self.assertEqual(
                ReactionEnsembleTest.reactant_coefficients[i],
                forward_reaction["reactant_coefficients"][i],
                msg="reactant coefficients not set correctly.")
        for i in range(len(forward_reaction["product_types"])):
            self.assertEqual(
                ReactionEnsembleTest.product_types[i],
                forward_reaction["product_types"][i],
                msg="product type not set correctly.")
        for i in range(len(forward_reaction["product_coefficients"])):
            self.assertEqual(
                ReactionEnsembleTest.product_coefficients[i],
                forward_reaction["product_coefficients"][i],
                msg="product coefficients not set correctly.")

        self.assertAlmostEqual(
            ReactionEnsembleTest.temperature,
            RE_status["temperature"],
            places=9,
            msg="reaction ensemble temperature not set correctly.")
        self.assertAlmostEqual(
            ReactionEnsembleTest.exclusion_radius,
            RE_status["exclusion_radius"],
            places=9,
            msg="reaction ensemble exclusion radius not set correctly.")

        self.assertAlmostEqual(
            ReactionEnsembleTest.volume,
            ReactionEnsembleTest.RE.get_volume(),
            places=9,
            msg="reaction ensemble volume not set correctly.")


if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
