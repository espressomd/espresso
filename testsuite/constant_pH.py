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
import sys
import os
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
    type_H = 5
    temperature = 1.0
    standard_pressure_in_simulation_units = 0.00108
    # avoid extreme regions in the titration curve e.g. via the choice
    # (np.random.random()-0.5)
    pKa_minus_pH = 1
    pH = 2  # or randomly via: 4*np.random.random()
    pKa = pKa_minus_pH + pH
    # could be in this test for example anywhere in the range 0.000001 ... 9
    K_HA_diss_apparent = 10**(-pKa)
    box_l = (N0 / c0)**(1.0 / 3.0)
    system = espressomd.System()
    system.seed = system.cell_system.get_state()['n_nodes'] * [2]
    np.random.seed(69) #make reaction code fully deterministic
    system.box_l = [box_l, box_l, box_l]
    system.cell_system.skin = 0.4
    system.time_step = 0.01
    RE = reaction_ensemble.ConstantpHEnsemble(
        standard_pressure=standard_pressure_in_simulation_units,
        temperature=1.0,
        exclusion_radius=1)

    @classmethod
    def setUpClass(cls):
        """Prepare a testsystem."""
        for i in range(0, 2 * cls.N0, 2):
            cls.system.part.add(id=i, pos=np.random.random(
                3) * cls.system.box_l, type=cls.type_A)
            cls.system.part.add(id=i + 1, pos=np.random.random(3) *
                                cls.system.box_l, type=cls.type_H)

        cls.RE.add(
            equilibrium_constant=cls.K_HA_diss_apparent, reactant_types=[
                cls.type_HA], reactant_coefficients=[1], product_types=[
                cls.type_A, cls.type_H], product_coefficients=[
                1, 1])
        cls.RE.set_default_charges(dictionary={cls.type_HA: 0, cls.type_A: -1, cls.type_H: +1})
        cls.RE.constant_pH = cls.pH

    @classmethod
    def ideal_degree_of_association(cls, pH):
        return 1 - 1.0 / (1 + 10**(cls.pKa - pH))

    def test_ideal_titration_curve(self):
        N0 = ReactionEnsembleTest.N0
        temperature = ReactionEnsembleTest.temperature
        type_A = ReactionEnsembleTest.type_A
        type_H = ReactionEnsembleTest.type_H
        type_HA = ReactionEnsembleTest.type_HA
        box_l = ReactionEnsembleTest.system.box_l
        system = ReactionEnsembleTest.system
        RE = ReactionEnsembleTest.RE
        """ chemical warmup in order to get to chemical equilibrium before starting to calculate the observable "degree of association" """
        for i in range(40 * N0):
            r = RE.reaction()

        volume = np.prod(self.system.box_l)  # cuboid box
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
        # note you cannot calculate the pH via -log10(<NH>/volume) in the
        # constant pH ensemble, since the volume is totally arbitrary and does
        # not influence the average number of protons
        pH = ReactionEnsembleTest.pH
        real_error_in_degree_of_association = abs(
            average_degree_of_association - ReactionEnsembleTest.ideal_degree_of_association(
                ReactionEnsembleTest.pH)) / ReactionEnsembleTest.ideal_degree_of_association(
            ReactionEnsembleTest.pH)
        print(average_degree_of_association)
        self.assertLess(
            real_error_in_degree_of_association,
            0.07,
            msg="Deviation to ideal titration curve for the given input parameters too large.")


if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
