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

"""Testmodule for the Widom Insertion.
"""
import os
import sys
import unittest as ut
import numpy as np
import espressomd  # pylint: disable=import-error
from espressomd import reaction_ensemble
from tests_common import lj_potential


@ut.skipIf(not espressomd.has_features(["LENNARD_JONES"]),
           "Features not available, skipping test!")
class WidomInsertionTest(ut.TestCase):

    """Test the implementation of the widom insertion.

       The excess chemical potential is calculated for identical particles in a 20 cubed box with a single particle, interacting via a LJ-potential (cut-off at 5 sigma)."""

    N0 = 1
    TEMPERATURE = 0.5
    TYPE_HA = 0
    CHARGE_HA = 0
    LJ_EPS = 1.0
    LJ_SIG = 1.0
    LJ_CUT = 5
    BOX_L = 2 * LJ_CUT
    LJ_SHIFT = lj_potential(LJ_CUT, LJ_EPS, LJ_SIG, LJ_CUT + 1, 0.0)

    radius = np.linspace(1e-10, LJ_CUT, 1000)
    # numerical integration for radii smaller than the cut-off in spherical
    # coordinates
    integrateUpToCutOff = 4 * np.pi * np.trapz(
        radius**2 * np.exp(-lj_potential(radius,
                                         LJ_EPS,
                                         LJ_SIG,
                                         LJ_CUT,
                                         LJ_SHIFT) / TEMPERATURE),
        x=radius)   
    # numerical solution for V_lj=0 => corresponds to the volume (as exp(0)=1)
    integreateRest = (BOX_L**3 - 4.0 / 3.0 * np.pi * LJ_CUT**3)
                      
    # calculate excess chemical potential of the system, see Frenkel Smith, p
    # 174. Note: He uses scaled coordinates, which is why we need to divide by
    # the Box volume
    target_mu_ex = -TEMPERATURE * \
        np.log((integrateUpToCutOff + integreateRest) / BOX_L**3)
                                       
    system = espressomd.System(box_l=np.ones(3) * BOX_L)
    system.cell_system.set_n_square()
    system.seed = system.cell_system.get_state()['n_nodes'] * [2]
    np.random.seed(69)  # make reaction code fully deterministic
    system.cell_system.skin = 0.4
    volume = np.prod(system.box_l)  # cuboid box
    
    Widom = reaction_ensemble.WidomInsertion(
        temperature=TEMPERATURE,
        exclusion_radius=np.nan)    

    def setUp(self):
        """Prepare a testsystem."""
        self.system.part.add(
            id=0,
            pos=0.5 *
            self.system.box_l,
            type=self.TYPE_HA)
        
        self.system.non_bonded_inter[self.TYPE_HA, self.TYPE_HA].lennard_jones.set_params(
            epsilon=self.LJ_EPS,
                                                           sigma=self.LJ_SIG,
                                                           cutoff=self.LJ_CUT, shift="auto")

        self.Widom.add_reaction(
            gamma=np.nan,
            reactant_types=[],
            reactant_coefficients=[],
            product_types=[self.TYPE_HA],
            product_coefficients=[1],
            default_charges={self.TYPE_HA: self.CHARGE_HA})

    def test_widom_insertion(self):   
        TYPE_HA = WidomInsertionTest.TYPE_HA
        system = WidomInsertionTest.system
        Widom = WidomInsertionTest.Widom
        target_mu_ex = WidomInsertionTest.target_mu_ex

        system.seed = system.cell_system.get_state()[
            'n_nodes'] * [np.random.randint(5)]
        num_samples = 100000
        for i in range(num_samples):
            Widom.measure_excess_chemical_potential(
                0)  # 0 for insertion reaction
        mu_ex = Widom.measure_excess_chemical_potential(0)
        deviation_mu_ex = abs(mu_ex[0] - target_mu_ex)

        # error
        self.assertLess(
            deviation_mu_ex - 1e-3,
            0.0,
            msg="\nExcess chemical potential for single LJ-particle computed via widom insertion gives a wrong value.\n"
            + "  average mu_ex: " + str(mu_ex[0])
            + "   mu_ex_std_err: " + str(mu_ex[1])
            + "  target_mu_ex: " + str(target_mu_ex)
        )

if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
