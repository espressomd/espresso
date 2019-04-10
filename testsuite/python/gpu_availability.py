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
# Tests particle property setters/getters
from __future__ import print_function
import unittest as ut
import espressomd


class GPUAvailability(ut.TestCase):

    """Tests consistency of gpu availability reporting."""

    def test(self):
        if "CUDA" in espressomd.features():
            system = espressomd.System(box_l=[1, 1, 1])
            self.assertEqual(
                system.cuda_init_handle.device_list != {},
                espressomd.gpu_available())
        else:
            self.assertFalse(espressomd.gpu_available())

if __name__ == "__main__":
    # print("Features: ", espressomd.features())
    ut.main()
