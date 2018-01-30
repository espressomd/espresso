#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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

    N0 = 40
    c0 = 0.00028
    type_HA = 0
    type_A = 1
    type_H = 2
    temperature = 1.0
    standard_pressure_in_simulation_units = 0.00108
    exclusion_radius = 1.0
    # could be in this test for example anywhere in the range 0.000001 ... 9,
    # chosen for example like 8.8*np.random.random()+0.2
    K_HA_diss = 8.8 * 0.5 + 0.2
    reactant_types = [type_HA]
    reactant_coefficients = [1]
    product_types = [type_A, type_H]
    product_coefficients = [1, 1]
    system = espressomd.System()
    system.seed = system.cell_system.get_state()['n_nodes'] * [2]
    np.random.seed(69) #make reaction code fully deterministic
    system.box_l = np.ones(3) * (N0 / c0)**(1.0 / 3.0)
    system.cell_system.skin = 0.4
    system.time_step = 0.01
    RE = reaction_ensemble.ReactionEnsemble(
        standard_pressure=standard_pressure_in_simulation_units,
        temperature=temperature,
        exclusion_radius=exclusion_radius)
    volume = np.prod(system.box_l)  # cuboid box

    @classmethod
    def setUpClass(cls):
        """Prepare a testsystem."""
        for i in range(0, 2 * cls.N0, 2):
            cls.system.part.add(id=i, pos=np.random.random(
                3) * cls.system.box_l, type=cls.type_A)
            cls.system.part.add(id=i + 1, pos=np.random.random(3) *
                                cls.system.box_l, type=cls.type_H)

        cls.RE.add(
            equilibrium_constant=cls.K_HA_diss,
            reactant_types=cls.reactant_types,
            reactant_coefficients=cls.reactant_coefficients,
            product_types=cls.product_types,
            product_coefficients=cls.product_coefficients)
        cls.RE.set_default_charges(dictionary={"0": 0, "1": -1, "2": +1})

    @classmethod
    def ideal_degree_of_association(cls, pK_a, pH):
        return 1 - 1.0 / (1 + 10**(pK_a - pH))

    def test_ideal_titration_curve(self):
        N0 = ReactionEnsembleTest.N0
        temperature = ReactionEnsembleTest.temperature
        type_A = ReactionEnsembleTest.type_A
        type_H = ReactionEnsembleTest.type_H
        type_HA = ReactionEnsembleTest.type_HA
        box_l = ReactionEnsembleTest.system.box_l
        standard_pressure_in_simulation_units = ReactionEnsembleTest.standard_pressure_in_simulation_units
        system = ReactionEnsembleTest.system
        K_HA_diss = ReactionEnsembleTest.K_HA_diss
        RE = ReactionEnsembleTest.RE
        """ chemical warmup in order to get to chemical equilibrium before starting to calculate the observable "degree of association" """
        RE.reaction(40 * N0)

        volume = ReactionEnsembleTest.volume
        average_NH = 0.0
        average_degree_of_association = 0.0
        num_samples = 1000
        for i in range(num_samples):
            RE.reaction()
            average_NH += system.number_of_particles(
                type=type_H)
            average_degree_of_association += system.number_of_particles(
                type=type_HA) / float(N0)
        average_NH /= num_samples
        average_degree_of_association /= num_samples
        pH = -np.log10(average_NH / volume)
        K_apparent_HA_diss = K_HA_diss * standard_pressure_in_simulation_units / temperature
        pK_a = -np.log10(K_apparent_HA_diss)
        print(average_degree_of_association)
        real_error_in_degree_of_association = abs(
            average_degree_of_association - ReactionEnsembleTest.ideal_degree_of_association(
                pK_a, pH)) / ReactionEnsembleTest.ideal_degree_of_association(
            pK_a, pH)
        self.assertLess(
            real_error_in_degree_of_association,
            0.07,
            msg="Deviation to ideal titration curve for the given input parameters too large.")

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

        self.assertAlmostEqual(
            ReactionEnsembleTest.standard_pressure_in_simulation_units,
            RE_status["standard_pressure"],
            places=9,
            msg="reaction ensemble standard_pressure not set correctly.")

if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
