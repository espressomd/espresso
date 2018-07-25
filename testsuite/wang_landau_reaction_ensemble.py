#
# Copyright (C) 2013,2014 The ESPResSo project
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
"""Testmodule for the Wang-Landau Reaction Ensemble.
"""
from __future__ import print_function
import numpy as np
import unittest as ut

import espressomd
from espressomd import code_info
from espressomd import analyze
from espressomd import integrate
from espressomd.interactions import *
from espressomd import reaction_ensemble
from espressomd import system


class ReactionEnsembleTest(ut.TestCase):

    """Test the core implementation of the wang_landau reaction ensemble.
    
    Create a harmonic bond between the two reacting particles. Therefore the
    potential energy is quadratic in the elongation of the bond and
    therefore the density of states is known as the one of the harmonic
    oscillator
    """

    # System parameters
    #
    box_l = 6 * np.sqrt(2)
    temperature = 1.0

    # Integration parameters
    #
    system = espressomd.System(box_l = [box_l, box_l, box_l])
    system.seed = system.cell_system.get_state()['n_nodes'] * [1234]
    np.random.seed(seed=system.seed)
    system.time_step = 0.01
    system.cell_system.skin = 0.4

    #
    # Setup System
    #

    N0 = 1  # number of titratable units
    K_diss = 0.0088

    system.part.add(id=0, pos=[0, 0, 0] * system.box_l, type=3)
    system.part.add(id=1, pos=[1.0, 1.0, 1.0] * system.box_l / 2.0, type=1)
    system.part.add(id=2, pos=np.random.random() * system.box_l, type=2)
    system.part.add(id=3, pos=np.random.random() * system.box_l, type=2)


    h = HarmonicBond(r_0=0, k=1)
    system.bonded_inter[0] = h
    system.part[0].add_bond((h, 1))
    RE = reaction_ensemble.WangLandauReactionEnsemble(temperature=temperature, exclusion_radius=0)
    RE.add_reaction(gamma=K_diss, reactant_types=[0], reactant_coefficients=[
           1], product_types=[1, 2], product_coefficients=[1, 1], default_charges={0: 0, 1: -1, 2: +1})
    system.setup_type_map([0, 1, 2, 3])
    # initialize wang_landau
    # generate preliminary_energy_run_results here, this should be done in a
    # seperate simulation without energy reweighting using the update energy
    # functions
    np.savetxt("energy_boundaries.dat", np.c_[
               [0, 1], [0, 0], [9, 9]],delimiter='\t' ,  header="nbar   E_potmin   E_potmax")

    RE.add_collective_variable_degree_of_association(
        associated_type=0, min=0, max=1, corresponding_acid_types=[0, 1])
    RE.add_collective_variable_potential_energy(
        filename="energy_boundaries.dat", delta=0.05)
    RE.set_wang_landau_parameters(
        final_wang_landau_parameter=1e-2, do_not_sample_reaction_partition_function=True, full_path_to_output_filename="WL_potential_out.dat")

    def test_wang_landau_output(self):
        while True:
            try:
                self.RE.reaction()
                for i in range(3):
                    self.RE.displacement_mc_move_for_particles_of_type(3)
            except reaction_ensemble.WangLandauHasConverged:  # only catch my exception
                break
        # test as soon as wang_landau has converged (throws exception then)
        nbars, Epots, WL_potentials = np.loadtxt(
            "WL_potential_out.dat", unpack=True)
        mask_nbar_0 = np.where(np.abs(nbars - 1.0) < 0.0001)
        Epots = Epots[mask_nbar_0]
        Epots = Epots[1:]
        WL_potentials = WL_potentials[mask_nbar_0]
        WL_potentials = WL_potentials[1:]

        expected_canonical_potential_energy = np.sum(np.exp(WL_potentials) * Epots * np.exp(
            -Epots / self.temperature)) / np.sum(np.exp(WL_potentials) * np.exp(-Epots / self.temperature))

        expected_canonical_squared_potential_energy = np.sum(np.exp(WL_potentials) * Epots**2 * np.exp(
            -Epots / self.temperature)) / np.sum(np.exp(WL_potentials) * np.exp(-Epots / self.temperature))

        expected_canonical_configurational_heat_capacity = expected_canonical_squared_potential_energy - \
            expected_canonical_potential_energy**2

        print(expected_canonical_potential_energy, expected_canonical_configurational_heat_capacity)

        # for the calculation regarding the analytical results which are
        # compared here, see Master Thesis Jonas Landsgesell p. 72
        self.assertAlmostEqual(
            expected_canonical_potential_energy - 1.5, 0.00, places=1,
                               msg="difference to analytical expected canonical potential energy too big")
        self.assertAlmostEqual(
            expected_canonical_configurational_heat_capacity - 1.5, 0.00, places=1,
                               msg="difference to analytical expected canonical configurational heat capacity too big")

if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
