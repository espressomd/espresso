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
"""Testmodule for the Wang-Landau Reaction Ensemble.
"""
import numpy as np
import unittest as ut

import espressomd
import espressomd.interactions
import espressomd.reaction_ensemble


class WangLandauReactionEnsembleTest(ut.TestCase):

    """Test the interface of the Wang-Landau reaction ensemble."""

    # System parameters
    #
    box_l = 6 * np.sqrt(2)
    temperature = 1.0

    # Integration parameters
    #
    system = espressomd.System(box_l=[box_l, box_l, box_l])
    np.random.seed(seed=42)
    system.time_step = 0.01
    system.cell_system.skin = 0
    system.cell_system.set_n_square(use_verlet_lists=False)

    #
    # Setup System
    #

    K_diss = 0.0088

    p1 = system.part.add(type=3, pos=[0, 0, 0])
    p2 = system.part.add(type=1, pos=system.box_l / 2.0)
    system.part.add(type=2, pos=np.random.random() * system.box_l)
    system.part.add(type=2, pos=np.random.random() * system.box_l)

    harmonic_bond = espressomd.interactions.HarmonicBond(r_0=0, k=1)
    system.bonded_inter[0] = harmonic_bond
    p1.add_bond((harmonic_bond, p2))

    WLRE = espressomd.reaction_ensemble.WangLandauReactionEnsemble(
        temperature=temperature, exclusion_radius=0, seed=86)
    WLRE.add_reaction(
        gamma=K_diss, reactant_types=[0], reactant_coefficients=[1],
        product_types=[1, 2], product_coefficients=[1, 1],
        default_charges={0: 0, 1: -1, 2: +1})
    system.setup_type_map([0, 1, 2, 3])
    # initialize wang_landau
    file_input = "energy_boundaries.dat"
    file_output = "WL_potential_out.dat"
    # generate preliminary_energy_run_results here, this should be done in a
    # separate simulation without energy reweighting using the update energy
    # functions
    np.savetxt(file_input, np.transpose([[0, 1], [0, 0], [9, 9]]),
               delimiter='\t', header="nbar   E_potmin   E_potmax")

    WLRE.add_collective_variable_degree_of_association(
        associated_type=0, min=0, max=1, corresponding_acid_types=[0, 1])
    WLRE.set_wang_landau_parameters(
        final_wang_landau_parameter=0.8 * 1e-2,
        do_not_sample_reaction_partition_function=True,
        full_path_to_output_filename=file_output)

    def test_01_energy_recording(self):
        self.WLRE.update_maximum_and_minimum_energies_at_current_state()
        self.WLRE.write_out_preliminary_energy_run_results()
        nbars, E_mins, E_maxs = np.loadtxt(
            "preliminary_energy_run_results", unpack=True)
        np.testing.assert_almost_equal(nbars, [0, 1])
        np.testing.assert_almost_equal(E_mins, [27.0, -10])
        np.testing.assert_almost_equal(E_maxs, [27.0, -10])

    def check_checkpoint(self, filename):
        # write first checkpoint
        self.WLRE.write_wang_landau_checkpoint()
        old_checkpoint = np.loadtxt(filename)

        # modify old_checkpoint in memory and in file (this destroys the
        # information contained in the checkpoint, but allows for testing of
        # the functions)
        modified_checkpoint = old_checkpoint
        modified_checkpoint[0] = 1
        np.savetxt(filename, modified_checkpoint)

        # check whether changes are carried out correctly
        self.WLRE.load_wang_landau_checkpoint()
        self.WLRE.write_wang_landau_checkpoint()
        new_checkpoint = np.loadtxt(filename)
        np.testing.assert_almost_equal(new_checkpoint, modified_checkpoint)

    def test_02_checkpointing(self):
        self.WLRE.add_collective_variable_potential_energy(
            filename=self.file_input, delta=0.05)

        # run MC for long enough to sample a non-trivial histogram
        for _ in range(1000):
            try:
                self.WLRE.reaction()
                for _ in range(2):
                    self.WLRE.displacement_mc_move_for_particles_of_type(3)
            except espressomd.reaction_ensemble.WangLandauHasConverged:
                break

        filenames = ["checkpoint_wang_landau_potential_checkpoint",
                     "checkpoint_wang_landau_histogram_checkpoint"]
        for filename in filenames:
            self.check_checkpoint(filename)


if __name__ == "__main__":
    ut.main()
