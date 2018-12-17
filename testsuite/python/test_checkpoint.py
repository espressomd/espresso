# Copyright (C) 2010-2018 The ESPResSo project
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
import subprocess
import unittest as ut
import numpy as np

import espressomd
import espressomd.checkpointing
import espressomd.virtual_sites
import tests_common

class CheckpointTest(ut.TestCase):

    @classmethod
    def setUpClass(self):
        checkpoint = espressomd.checkpointing.Checkpoint(
            checkpoint_id="mycheckpoint",
            checkpoint_path="@CMAKE_CURRENT_BINARY_DIR@")
        checkpoint.load(0)
        if espressomd.has_features('LB'):
            self.lbf = system.actors[0]
            self.lbf.load_checkpoint("@CMAKE_CURRENT_BINARY_DIR@/lb.cpt", 1)

    @ut.skipIf(not espressomd.has_features('LB'),
               "Skipping test due to missing features.")
    def test_LB(self):
        np.testing.assert_almost_equal(
            np.copy(self.lbf[1, 1, 1].velocity), np.array([0.1, 0.2, 0.3]))
        state = self.lbf.get_params()
        reference = {'agrid': 0.5, 'visc': 1.3,
                     'dens': 1.5, 'tau': 0.01, 'fric': 2.0}
        self.assertDictContainsSubset(reference, state)

    def test_variables(self):
        self.assertEqual(system.cell_system.skin, 0.1)
        self.assertEqual(system.time_step, 0.01)
        self.assertEqual(system.min_global_cut, 2.0)

    def test_part(self):
        np.testing.assert_allclose(
            np.copy(system.part[0].pos), np.array([1.0, 2.0, 3.0]))
        np.testing.assert_allclose(
            np.copy(system.part[1].pos), np.array([1.0, 1.0, 2.0]))
        np.testing.assert_allclose(np.copy(system.part[0].f), particle_force0)
        np.testing.assert_allclose(np.copy(system.part[1].f), particle_force1)

    def test_thermostat(self):
        self.assertEqual(system.thermostat.get_state()[0]['type'], 'LANGEVIN')
        self.assertEqual(system.thermostat.get_state()[0]['kT'], 1.0)
        np.testing.assert_array_equal(system.thermostat.get_state()[
            0]['gamma'], np.array([2.0, 2.0, 2.0]))

    @ut.skipIf(
        not espressomd.has_features(
            ['LENNARD_JONES']),
        "Cannot test for Lennard-Jones checkpointing because feature not compiled in.")
    def test_non_bonded_inter(self):
        state = system.non_bonded_inter[
            0, 0].lennard_jones._get_params_from_es_core()
        state2 = system.non_bonded_inter[
            3, 0].lennard_jones._get_params_from_es_core()
        reference = {'shift': 0.1, 'sigma': 1.3, 'epsilon':
                     1.2, 'cutoff': 2.0, 'offset': 0.0, 'min': 0.0}
        reference2 = {'shift': 0.1, 'sigma': 1.7, 'epsilon':
                      1.2, 'cutoff': 2.0, 'offset': 0.0, 'min': 0.0}
        self.assertEqual(
            len(set(state.items()) & set(reference.items())), len(reference))
        self.assertEqual(len(set(state2.items()) & set(
            reference2.items())), len(reference2))

    def test_bonded_inter(self):
        state = system.part[1].bonds[0][0].params
        reference = {'r_0': 0.0, 'k': 1.0}
        self.assertEqual(
            len(set(state.items()) & set(reference.items())), len(reference))

    @ut.skipIf(
        not espressomd.has_features(
            ['VIRTUAL_SITES', 'VIRTUAL_SITES_RELATIVE']),
        "Cannot test for virtual site checkpointing because feature not compiled in.")
    def test_virtual_sites(self):
        self.assertEqual(system.part[1].virtual, 1)
        self.assertTrue(
            isinstance(
                system.virtual_sites,
                espressomd.virtual_sites.VirtualSitesRelative))

    def test_mean_variance_calculator(self):
        np.testing.assert_array_equal(
            acc.get_mean(), np.array([1.0, 1.5, 2.0, 1.0, 1.0, 2.0]))
        np.testing.assert_array_equal(
            acc.get_variance(), np.array([0., 0.5, 2., 0., 0., 0.]))

    @ut.skipIf(
        not espressomd.has_features(
            ['ELECTROSTATICS']),
        "Cannot test for P3M checkpointing because feature not compiled in.")
    def test_p3m(self):
        self.assertTrue(any(isinstance(actor, espressomd.electrostatics.P3M)
                            for actor in system.actors.active_actors))
    
    @ut.skipIf(not espressomd.has_features("COLLISION_DETECTION"), "skipped for missing features")
    def test_collision_detection(self):
        coldet = system.collision_detection
        self.assertEqual(coldet.mode, "bind_centers")
        self.assertAlmostEqual(coldet.distance, 0.11, delta=1E-9)
        self.assertTrue(coldet.bond_centers, system.bonded_inter[0])

    @ut.skipIf(not espressomd.has_features("EXCLUSIONS"), "Skipped because feature EXCLUSIONS missing.")
    def test_exclusions(self):
        self.assertTrue(tests_common.lists_contain_same_elements(system.part[0].exclusions, [2]))
        self.assertTrue(tests_common.lists_contain_same_elements(system.part[1].exclusions, [2]))
        self.assertTrue(tests_common.lists_contain_same_elements(system.part[2].exclusions, [0,1]))


if __name__ == '__main__':
    ut.main()
