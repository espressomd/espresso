import subprocess
import unittest as ut
import numpy as np

import espressomd
import espressomd.checkpointing
import espressomd.virtual_sites



class CheckpointTest(ut.TestCase):
    @classmethod
    def setUpClass(self):
        # Write checkpoint.
        p = subprocess.Popen(['@CMAKE_BINARY_DIR@/pypresso', '@CMAKE_CURRENT_BINARY_DIR@/save_checkpoint.py'])
        p.wait()
        system = espressomd.System(box_l=[10.0, 10.0, 10.0])
        checkpoint = espressomd.checkpointing.Checkpointing(checkpoint_id="mycheckpoint", checkpoint_path="@CMAKE_CURRENT_BINARY_DIR@")
        checkpoint.load(0)

    def test_variables(self):
        self.assertEqual(system.cell_system.skin, 0.4)
        self.assertEqual(system.time_step, 0.01)
        self.assertEqual(system.min_global_cut, 2.0)

    def test_part(self):
        np.testing.assert_array_equal(np.copy(system.part[0].pos), np.array([1.0, 1.0, 1.0]))
        np.testing.assert_array_equal(np.copy(system.part[1].pos), np.array([1.0, 1.0, 2.0]))

    def test_thermostat(self):
        self.assertEqual(system.thermostat.get_state()[0]['type'], 'LANGEVIN')
        self.assertEqual(system.thermostat.get_state()[0]['kT'], 1.0)
        np.testing.assert_array_equal(system.thermostat.get_state()[0]['gamma'], np.array([2.0, 2.0, 2.0]))

    def test_non_bonded_inter(self):
        state = system.non_bonded_inter.__getstate__()[0][0]
        reference = {'shift': 0.1, 'sigma': 1.3, 'epsilon': 1.2, 'cutoff': 2.0, 'type_name': 'LennardJones', 'offset': 0.0, 'min': 0.0}
        self.assertEqual(len(set(state.items()) & set(reference.items())), 7)

    @ut.skipIf(not espressomd.has_features(['VIRTUAL_SITES', 'VIRTUAL_SITES_RELATIVE']),
               "Cannot test for virtual site checkpointing because feature not compiled in.")
    def test_virtual_sites(self):
        self.assertEqual(system.part[1].virtual, 1)
        self.assertTrue(isinstance(system.virtual_sites, espressomd.virtual_sites.VirtualSitesRelative))

if __name__ == '__main__':
    ut.main()
