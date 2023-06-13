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
import importlib_wrapper
import numpy as np

np.random.seed(42)

tutorial, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@TUTORIALS_DIR@/raspberry_electrophoresis/raspberry_electrophoresis.py",
    box_l=16., num_iterations=100, num_steps_per_iteration=80)


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    '''Check the raspberry travels a longer distance on the x-axis'''
    system = tutorial.system

    def test_steepest_descent_convergence(self):
        self.assertLess(tutorial.force_max, 10.)

    def test_trajectory_sample(self):
        trajectory = np.loadtxt('posVsTime_sample.dat')[:, 1:4]
        x, y, z = np.abs(trajectory[-1, :] - trajectory[0, :])
        self.assertGreater(x, y)
        self.assertGreater(x, z)

    def test_trajectory_simulated(self):
        trajectory = np.loadtxt('posVsTime.dat')[:, 1:4]
        x, y, z = np.abs(trajectory[-1, :] - trajectory[0, :])
        self.assertGreater(x, y)
        self.assertGreater(x, z)


if __name__ == "__main__":
    ut.main()
