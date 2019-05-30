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
import sys
import subprocess
import unittest as ut
import numpy as np

import espressomd
import espressomd.checkpointing
import espressomd.virtual_sites
import tests_common

modes = {x for mode in set("@TEST_COMBINATION@".upper().split('-'))
         for x in [mode, mode.split('.')[0]]}

LB = ('LB.CPU' in modes or
      espressomd.gpu_available() and espressomd.has_features('CUDA') and 'LB.GPU' in modes)

EK = (espressomd.gpu_available() and espressomd.has_features('ELECTROKINETICS') and 'EK.GPU'
      in modes)


class CheckpointTest(ut.TestCase):

    @classmethod
    def setUpClass(self):
        self.checkpoint = espressomd.checkpointing.Checkpoint(
            checkpoint_id="mycheckpoint_@TEST_COMBINATION@_@TEST_BINARY@".replace(
                '.', '__'),
            checkpoint_path="@CMAKE_CURRENT_BINARY_DIR@")
        self.checkpoint.load(0)

    @ut.skipIf(not LB, "Skipping test due to missing features.")
    def test_LB(self):
        lbf = system.actors[0]
        cpt_mode = int("@TEST_BINARY@")
        cpt_path = self.checkpoint.checkpoint_dir + "/lb{}.cpt"
        with self.assertRaises(RuntimeError):
            lbf.load_checkpoint(cpt_path.format("-corrupted"), cpt_mode)
        assertRaisesRegex = self.assertRaisesRegexp if sys.version_info < (
            3, 2) else self.assertRaisesRegex
        with assertRaisesRegex(RuntimeError, 'grid dimensions mismatch'):
            lbf.load_checkpoint(cpt_path.format("-wrong-boxdim"), cpt_mode)
        lbf.load_checkpoint(cpt_path.format(""), cpt_mode)
        precision = 9 if "LB.CPU" in modes else 5
        m = np.pi / 12
        nx = int(np.round(system.box_l[0] / lbf.get_params()["agrid"]))
        ny = int(np.round(system.box_l[1] / lbf.get_params()["agrid"]))
        nz = int(np.round(system.box_l[2] / lbf.get_params()["agrid"]))
        grid_3D = np.fromfunction(
            lambda i, j, k: np.cos(i * m) * np.cos(j * m) * np.cos(k * m),
            (nx, ny, nz), dtype=float)
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    np.testing.assert_almost_equal(
                        np.copy(lbf[i, j, k].population),
                                grid_3D[i, j, k] * np.arange(1, 20),
                                decimal=precision)
        state = lbf.get_params()
        reference = {'agrid': 0.5, 'visc': 1.3,
                     'dens': 1.5, 'tau': 0.01}
        for key, val in reference.items():
            self.assertTrue(key in state)
            self.assertAlmostEqual(reference[key], state[key], delta=1E-9)

    @ut.skipIf(not EK, "Skipping test due to missing features.")
    def test_EK(self):
        ek = system.actors[0]
        ek_species = ek.get_params()['species'][0]
        cpt_path = self.checkpoint.checkpoint_dir + "/ek"
        ek.load_checkpoint(cpt_path)
        precision = 5
        m = np.pi / 12
        nx = int(np.round(system.box_l[0] / ek.get_params()["agrid"]))
        ny = int(np.round(system.box_l[1] / ek.get_params()["agrid"]))
        nz = int(np.round(system.box_l[2] / ek.get_params()["agrid"]))
        grid_3D = np.fromfunction(
            lambda i, j, k: np.cos(i * m) * np.cos(j * m) * np.cos(k * m),
            (nx, ny, nz), dtype=float)
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    np.testing.assert_almost_equal(
                        np.copy(ek_species[i, j, k].density),
                                grid_3D[i, j, k],
                                decimal=precision)
        state = ek.get_params()
        reference = {'agrid': 0.5, 'lb_density': 26.15,
                     'viscosity': 1.7, 'friction': 0.0,
                     'T': 1.1, 'prefactor': 0.88, 'stencil': "linkcentered"}
        for key, val in reference.items():
            self.assertTrue(key in state)
            self.assertAlmostEqual(reference[key], state[key], delta=1E-5)
        state_species = ek_species.get_params()
        reference_species = {'density': 0.4, 'D': 0.02, 'valency': 0.3}
        for key, val in reference_species.items():
            self.assertTrue(key in state_species)
            self.assertAlmostEqual(
                reference_species[key],
                state_species[key],
                delta=1E-5)
        self.assertAlmostEqual(
            state_species['ext_force_density'][0],
            0.01,
            delta=1E-5)
        self.assertAlmostEqual(
            state_species['ext_force_density'][1],
            -0.08,
            delta=1E-5)
        self.assertAlmostEqual(
            state_species['ext_force_density'][2],
            0.06,
            delta=1E-5)

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
        self.assertEqual(system.thermostat.get_state()[0]['seed'], 42)
        np.testing.assert_array_equal(system.thermostat.get_state()[
            0]['gamma'], np.array([2.0, 2.0, 2.0]))

    @ut.skipIf(not espressomd.has_features('LENNARD_JONES'),
               "Skipping test due to missing features.")
    @ut.skipIf('LJ' not in modes,
               "Skipping test due to missing combination.")
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

    @ut.skipIf(not espressomd.has_features(['VIRTUAL_SITES',
                                            'VIRTUAL_SITES_RELATIVE']),
               "Skipping test due to missing features.")
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

    @ut.skipIf(not espressomd.has_features('ELECTROSTATICS'),
               "Skipping test due to missing features.")
    @ut.skipIf('P3M.CPU' not in modes,
               "Skipping test due to missing combination.")
    def test_p3m(self):
        self.assertTrue(any(isinstance(actor, espressomd.electrostatics.P3M)
                            for actor in system.actors.active_actors))
    
    @ut.skipIf(not espressomd.has_features('COLLISION_DETECTION'),
               "Skipping test due to missing features.")
    def test_collision_detection(self):
        coldet = system.collision_detection
        self.assertEqual(coldet.mode, "bind_centers")
        self.assertAlmostEqual(coldet.distance, 0.11, delta=1E-9)
        self.assertTrue(coldet.bond_centers, system.bonded_inter[0])

    @ut.skipIf(not espressomd.has_features('EXCLUSIONS'),
               "Skipping test due to missing features.")
    def test_exclusions(self):
        self.assertTrue(tests_common.lists_contain_same_elements(
            system.part[0].exclusions, [2]))
        self.assertTrue(tests_common.lists_contain_same_elements(
            system.part[1].exclusions, [2]))
        self.assertTrue(tests_common.lists_contain_same_elements(
            system.part[2].exclusions, [0, 1]))


if __name__ == '__main__':
    ut.main()
