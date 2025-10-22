#
# Copyright (C) 2024 The ESPResSo project
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
import espressomd


class Test(ut.TestCase):

    def test_system(self):
        with self.assertRaisesRegex(ValueError, "Required argument 'box_l' not provided"):
            espressomd.System()
        with self.assertRaisesRegex(ValueError, "Property 'unknown' cannot be set via argument to System class"):
            espressomd.System(box_l=[1., 1., 1.], unknown=1)
        system = espressomd.System(box_l=[1., 1., 1.], min_global_cut=0.01)
        self.assertEqual(system.min_global_cut, 0.01)
        self.assertEqual(system.time_step, -1.)
        with self.assertRaisesRegex(RuntimeError, "You can only have one instance of the system class at a time"):
            espressomd.System(box_l=[1., 1., 1.])


if __name__ == "__main__":
    ut.main()
