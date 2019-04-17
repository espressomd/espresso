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
class ReactionEnsembleTest(ut.TestCase):

    """Test the implementation of the widom insertion (part of the reaction ensemble)."""

    # The excess chemical potential is calculated for identical particles into
    # a 20 cubed box with a single particle, interacting via a LJ-potential
    # (cut-off at 5 sigma).
    N0 = 1
    temperature = 0.5
    type_HA = 0
    charge_HA = 0
    lj_eps = 1.0
    lj_sig = 1.0
    lj_cut = 5.0
    box_l = 20.0
    lj_shift = lj_potential(lj_cut, lj_eps, lj_sig, lj_cut + 1, 0.0)

    def lj_potential_arr(r, lj_eps, lj_sig, lj_shift):
        """lennard jones potential (that works for array input)"""
        r6 = r**(-6)
        r12 = r6**2
        ljsig6 = lj_sig**6
        ljsig12 = ljsig6**2
        return 4.0 * lj_eps * (ljsig12 * r12 - ljsig6 * r6) + lj_shift

    r = np.linspace(1e-10, lj_cut, 1000)
    integrateSphere = 4 * np.pi * np.trapz(
        r**2 * np.exp(-lj_potential_arr(r,
                                        lj_eps,
                                        lj_sig,
                                        lj_shift) / temperature),
        x=r)  # numerical solution 
    integreateRest = (box_l**3 - 4.0 / 3.0 * np.pi * lj_cut**3)
                      # numerical solution for V_lj=0 => corresponds to the
                      # volume (as exp(0)=1)

    target_mu_ex = -temperature * \
        np.log((integrateSphere + integreateRest) / box_l**3)
                                       # calculate excess chemical potential of
                                       # the system from the integral values

    system = espressomd.System(box_l=np.ones(3) * box_l)
    system.seed = system.cell_system.get_state()['n_nodes'] * [2]
    np.random.seed(69)  # make reaction code fully deterministic
    system.cell_system.skin = 0.4
    volume = np.prod(system.box_l)  # cuboid box
    
    unimportant_exclusion_radius = np.nan
    Widom = reaction_ensemble.WidomInsertion(
        temperature=temperature,
        exclusion_radius=unimportant_exclusion_radius)    

    @classmethod
    def setUpClass(cls):
        """Prepare a testsystem."""
        cls.system.part.add(id=0, pos=0.5 * cls.system.box_l, type=cls.type_HA)
        
        cls.system.non_bonded_inter[cls.type_HA, cls.type_HA].lennard_jones.set_params(
            epsilon=cls.lj_eps,
                                                           sigma=cls.lj_sig,
                                                           cutoff=cls.lj_cut, shift="auto")

        unimportant_K_diss = 0.8888
        cls.Widom.add_reaction(
            gamma=unimportant_K_diss,
            reactant_types=[],
            reactant_coefficients=[],
            product_types=[cls.type_HA],
            product_coefficients=[1],
            default_charges={cls.type_HA: cls.charge_HA})

    def test_widom_insertion(self):   
        type_HA = ReactionEnsembleTest.type_HA
        system = ReactionEnsembleTest.system
        Widom = ReactionEnsembleTest.Widom
        target_mu_ex = ReactionEnsembleTest.target_mu_ex

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
            deviation_mu_ex - 0.00015,
            0.0,
            msg="\nExcess chemical potential for single LJ-particle computed via widom insertion gives a wrong value.\n"
            + "  average mu_ex: " + str(mu_ex[0])
            + "   mu_ex_std_err: " + str(mu_ex[1])
            + "  target_mu_ex: " + str(target_mu_ex)
        )

if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
