# Copyright (C) 2019 The ESPResSo project
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

tutorial, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@TUTORIALS_DIR@/05-raspberry_electrophoresis/05-raspberry_electrophoresis.py",
    gpu=True, box_l=20., E=0.25, num_iterations=200, num_steps_per_iteration=250)


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    system = tutorial.system

    def test_trajectory_sample(self):
        trajectory = np.loadtxt('posVsTime_sample.dat')[:, 1:4]
        # the raspberry should have traveled mostly on the x-axis
        dist = np.abs(trajectory[-1, :] - trajectory[0, :])
        self.assertGreater(dist[0], dist[1])
        self.assertGreater(dist[0], dist[2])

    def test_trajectory_simulated(self):
        trajectory = np.loadtxt('posVsTime.dat')[:, 1:4]
        # the raspberry should have traveled mostly on the x-axis,
        # but due to insufficient sampling, it's not always the case
        dist = np.abs(trajectory[-1, :] - trajectory[0, :])
        self.assertGreater(dist[0], np.min(dist[1:]))


if __name__ == "__main__":
    ut.main()
