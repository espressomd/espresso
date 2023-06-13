# Copyright (C) 2019-2022 The ESPResSo project
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

import unittest as ut
import numpy as np
import importlib_wrapper


sample, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "/home/thilo/code/espresso2/espresso/testsuite/scripts/samples/local_samples/load_checkpoint.py")


@skipIfMissingFeatures
class Sample(ut.TestCase):
    system = sample.system

    def test_file_generation(self):
        self.assertEqual(set(sample.checkpoint.get_registered_objects()),
                         {'myvar', 'system', 'p3m'})
        self.assertEqual(sample.myvar, "some script variable (updated value)")

    def test_trajectory_reproducibility(self):
        self.assertTrue(sample.p3m.is_tuned)
        np.testing.assert_array_less(sample.forces_diff, 1e-16)


if __name__ == "__main__":
    ut.main()
