import subprocess
import unittest as ut
import numpy as np

import espressomd
import espressomd.checkpointing
import espressomd.virtual_sites



class CheckpointTest(ut.TestCase):
    @classmethod
    def setUpClass(self):
        checkpoint = espressomd.checkpointing.Checkpoint(checkpoint_id="mycheckpoint", checkpoint_path="@CMAKE_CURRENT_BINARY_DIR@")
        checkpoint.load(0)

    def test_variables(self):
        self.assertEqual(system.cell_system.skin, 0.4)
        self.assertEqual(system.time_step, 0.01)
        self.assertEqual(system.min_global_cut, 2.0)

    def test_part(self):
        np.testing.assert_allclose(np.copy(system.part[0].pos), np.array([1.0, 2.0, 3.0]))
        np.testing.assert_allclose(np.copy(system.part[1].pos), np.array([1.0, 1.0, 2.0]))
        np.testing.assert_allclose(np.copy(system.part[0].f), particle_force0)
        np.testing.assert_allclose(np.copy(system.part[1].f), particle_force1)

    def test_thermostat(self):
        self.assertEqual(system.thermostat.get_state()[0]['type'], 'LANGEVIN')
        self.assertEqual(system.thermostat.get_state()[0]['kT'], 1.0)
        np.testing.assert_array_equal(system.thermostat.get_state()[0]['gamma'], np.array([2.0, 2.0, 2.0]))

    @ut.skipIf(not espressomd.has_features(['LENNARD_JONES']),
               "Cannot test for Lennard-Jones checkpointing because feature not compiled in.")
    def test_non_bonded_inter(self):
        state = system.non_bonded_inter[0, 0].lennard_jones._get_params_from_es_core()
        state2 = system.non_bonded_inter[3, 0].lennard_jones._get_params_from_es_core()
        reference = {'shift': 0.1, 'sigma': 1.3, 'epsilon': 1.2, 'cutoff': 2.0, 'offset': 0.0, 'min': 0.0}
        reference2 = {'shift': 0.1, 'sigma': 1.7, 'epsilon': 1.2, 'cutoff': 2.0, 'offset': 0.0, 'min': 0.0}
        self.assertEqual(len(set(state.items()) & set(reference.items())), len(reference))
        self.assertEqual(len(set(state2.items()) & set(reference2.items())), len(reference2))

    def test_bonded_inter(self):
        state = system.part[1].bonds[0][0].params
        reference = {'r_0': 0.0, 'k': 1.0}
        self.assertEqual(len(set(state.items()) & set(reference.items())), len(reference))

    @ut.skipIf(not espressomd.has_features(['VIRTUAL_SITES', 'VIRTUAL_SITES_RELATIVE']),
               "Cannot test for virtual site checkpointing because feature not compiled in.")
    def test_virtual_sites(self):
        self.assertEqual(system.part[1].virtual, 1)
        self.assertTrue(isinstance(system.virtual_sites, espressomd.virtual_sites.VirtualSitesRelative))

    def test_mean_variance_calculator(self):
        np.testing.assert_array_equal(acc.get_mean(), np.array([1.0, 1.5, 2.0, 1.0, 1.0, 2.0]))
        np.testing.assert_array_equal(acc.get_variance(), np.array([0. ,  0.5,  2. ,  0. ,  0. ,  0.]))
    
    @ut.skipIf(not espressomd.has_features(['ELECTROSTATICS']),
              "Cannot test for P3M checkpointing because feature not compiled in.")  
    def test_p3m(self):
        self.assertTrue(isinstance(system.actors.active_actors[0], espressomd.electrostatics.P3M))

if __name__ == '__main__':
    ut.main()
