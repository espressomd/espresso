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
import unittest as ut
import unittest_decorators as utx
import espressomd


class GPUAvailability(ut.TestCase):

    """Tests consistency of GPU availability reporting."""
    system = espressomd.System(box_l=[1, 1, 1])

    def test(self):
        if espressomd.has_features("CUDA"):
            self.assertEqual(self.system.cuda_init_handle.list_devices() != {},
                             espressomd.gpu_available())
            self.assertEqual(
                self.system.cuda_init_handle.list_devices_properties() != {},
                espressomd.gpu_available())
        else:
            self.assertFalse(espressomd.gpu_available())

    @utx.skipIfMissingFeatures("CUDA")
    def test_exceptions(self):
        error_msg = 'CUDA error: '
        if espressomd.gpu_available():
            n_gpus = len(self.system.cuda_init_handle.list_devices())
            with self.assertRaisesRegex(RuntimeError, error_msg):
                self.system.cuda_init_handle.device = n_gpus + 1
        else:
            with self.assertRaisesRegex(RuntimeError, error_msg):
                self.system.cuda_init_handle.device
            with self.assertRaisesRegex(RuntimeError, error_msg):
                self.system.cuda_init_handle.device = 0

    @utx.skipIfMissingGPU()
    def test_list_devices(self):
        # check if GPU properties can be queried
        device_list = self.system.cuda_init_handle.list_devices()
        device_list_p = self.system.cuda_init_handle.list_devices_properties()
        self.assertEqual(len(device_list_p), 1)
        device_list_p_head = list(device_list_p.values())[0]
        dev_keys = {'name', 'compute_capability', 'cores', 'total_memory'}
        # check both dicts agree
        self.assertEqual(device_list.keys(), device_list_p_head.keys())
        for dev_id in device_list:
            self.assertEqual(device_list_p_head[dev_id].keys(), dev_keys)
            self.assertEqual(
                device_list_p_head[dev_id]['name'],
                device_list[dev_id])
        # check the currently active GPU
        dev_id = self.system.cuda_init_handle.device
        self.assertIn(dev_id, device_list_p_head)
        device = device_list_p_head[dev_id]
        self.assertGreater(device['cores'], 0)
        self.assertGreater(device['total_memory'], 0)
        self.assertGreaterEqual(device['compute_capability'][0], 3)
        self.assertGreaterEqual(device['compute_capability'][1], 0)


if __name__ == "__main__":
    ut.main()
