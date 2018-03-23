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
        checkpoint = espressomd.checkpointing.Checkpointing(checkpoint_id="mycheckpoint", checkpoint_path="@CMAKE_CURRENT_BINARY_DIR@")
        checkpoint.load(0)
        system = espressomd.System(box_l=[10.0, 10.0, 10.0])
        system.cell_system.skin = skin
        system.time_step = time_step

    def test_variables(self):
        self.assertEqual(skin, 0.4)
        self.assertEqual(time_step, 0.01)

    def test_part(self):
        np.testing.assert_array_equal(np.copy(system.part[0].pos), np.array([1.0, 1.0, 1.0]))
        np.testing.assert_array_equal(np.copy(system.part[1].pos), np.array([1.0, 1.0, 2.0]))

    @ut.skipIf(not espressomd.has_features(['VIRTUAL_SITES', 'VIRTUAL_SITES_RELATIVE']),
               "Cannot test for virtual site checkpointing because feature not compiled in.")
    def test_virtual_sites(self):
        self.assertEqual(system.part[1].virtual, 1)
        self.assertTrue(isinstance(system.virtual_sites, espressomd.virtual_sites.VirtualSitesRelative))

if __name__ == '__main__':
    ut.main()
