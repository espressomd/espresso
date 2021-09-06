#
# Copyright (C) 2013-2019 The ESPResSo project
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
import unittest as ut
import numpy as np
import espressomd
import espressomd.reaction_ensemble


class ReactionEnsembleTest(ut.TestCase):

    """Test the core implementation of the reaction ensemble."""

    # The reaction ensemble follows the ideal titration curve only if N>>1,
    # Ideal curve is derived in the grandcanonical ensemble and for low N
    # there are systematic deviations caused by differences between the
    # ensembles. This is not an error but a fundamental difference (various
    # ensembles are equivalent only in the thermodynamic limit N \to \infty)
    N0 = 40
    c0 = 0.00028
    types = {
        "HA": 0,
        "A-": 1,
        "H+": 2,
    }
    charge_dict = {
        0: 0,
        1: -1,
        2: +1,
    }
    target_alpha = 0.6
    # We get best statistics at alpha=0.5 Then the test is less sensitive to
    # the exact sequence of random numbers and does not require hard-coded
    # output values
    temperature = 1.0
    exclusion_radius = 1.0
    # could be in this test for example anywhere in the range 0.000001 ... 9,
    reactant_types = [types["HA"]]
    reactant_coefficients = [1]
    product_types = [types["A-"], types["H+"]]
    product_coefficients = [1, 1]
    nubar = 1
    system = espressomd.System(box_l=np.ones(3) * (N0 / c0)**(1.0 / 3.0))
    np.random.seed(69)  # make reaction code fully deterministic
    system.cell_system.skin = 0.4
    volume = system.volume()
    # Calculate gamma which should lead to target_alpha with given N0 and V
    # Requires N0>>1, otherwise discrete character of N changes the statistics (N>20 should suffice)
    # gamma = prod_i (N_i / V) = alpha^2 N0 / (1-alpha)*V**(-nubar)
    # degree of dissociation alpha = N_A / N_HA = N_H / N_0
    gamma = target_alpha**2 / (1. - target_alpha) * N0 / (volume**nubar)
    RE = espressomd.reaction_ensemble.ReactionEnsemble(
        kT=temperature,
        exclusion_radius=exclusion_radius, seed=12)

    @classmethod
    def setUpClass(cls):
        cls.system.part.add(
            pos=np.random.random((2 * cls.N0, 3)) * cls.system.box_l,
            type=cls.N0 * [cls.types["A-"], cls.types["H+"]])

        cls.RE.add_reaction(
            gamma=cls.gamma,
            reactant_types=cls.reactant_types,
            reactant_coefficients=cls.reactant_coefficients,
            product_types=cls.product_types,
            product_coefficients=cls.product_coefficients,
            default_charges=cls.charge_dict,
            check_for_electroneutrality=True)

    def test_ideal_titration_curve(self):
        N0 = ReactionEnsembleTest.N0
        types = ReactionEnsembleTest.types
        system = ReactionEnsembleTest.system
        gamma = ReactionEnsembleTest.gamma

        RE = ReactionEnsembleTest.RE
        target_alpha = ReactionEnsembleTest.target_alpha

        # Set the hidden particle type to the lowest possible number to speed
        # up the simulation
        RE.set_non_interacting_type(max(types.values()) + 1)

        # chemical warmup - get close to chemical equilibrium before we start
        # sampling
        RE.reaction(5 * N0)

        average_NH = 0.0
        average_NHA = 0.0
        average_NA = 0.0
        num_samples = 300
        for _ in range(num_samples):
            RE.reaction(10)
            average_NH += system.number_of_particles(type=types["H+"])
            average_NHA += system.number_of_particles(type=types["HA"])
            average_NA += system.number_of_particles(type=types["A-"])
        average_NH /= num_samples
        average_NA /= num_samples
        average_NHA /= num_samples
        average_alpha = average_NA / float(N0)
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
            + f"  gamma: {gamma:.3e}"
            + f"  average NH: {average_NH:.1f}"
            + f"  average NA: {average_NA:.1f}"
            + f"  average NHA: {average_NHA:.1f}"
            + f"  average alpha: {average_alpha:.3f}"
            + f"  target alpha: {target_alpha:.3f}"
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
            RE_status["kT"],
            places=9,
            msg="reaction ensemble kT not set correctly.")
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

    def test_change_reaction_constant(self):
        RE = ReactionEnsembleTest.RE
        new_reaction_constant = 634.0
        RE.change_reaction_constant(0, new_reaction_constant)
        RE_status = RE.get_status()
        forward_reaction = RE_status["reactions"][0]
        backward_reaction = RE_status["reactions"][1]
        self.assertEqual(
            new_reaction_constant,
            forward_reaction["gamma"],
            msg="new reaction constant was not set correctly.")
        self.assertEqual(
            1.0 / new_reaction_constant,
            backward_reaction["gamma"],
            msg="new reaction constant was not set correctly.")
        RE.change_reaction_constant(0, ReactionEnsembleTest.gamma)

    def test_delete_reaction(self):
        RE = ReactionEnsembleTest.RE
        RE.add_reaction(
            gamma=1,
            reactant_types=[5],
            reactant_coefficients=[1],
            product_types=[2, 3, 4],
            product_coefficients=[1, 4, 3],
            default_charges={5: 0, 2: 0, 3: 0, 4: 0},
            check_for_electroneutrality=True)
        nr_reactions_after_addition = len(RE.get_status()["reactions"])
        RE.delete_reaction(1)
        nr_reactions_after_deletion = len(RE.get_status()["reactions"])
        self.assertEqual(
            2,
            nr_reactions_after_addition - nr_reactions_after_deletion,
            msg="the difference in single reactions does not match,\
            deleting a full reaction (back and forward direction)\
            should result in deleting two single reactions.")


if __name__ == "__main__":
    ut.main()
