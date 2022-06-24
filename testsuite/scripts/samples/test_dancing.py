# Copyright (C) 2020-2022 The ESPResSo project
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

sd_method = "--fts" if "fts" in "@TEST_SUFFIX@" else "--ft"
if sd_method == "--fts":
    y_min = -555
    intsteps = 5500
else:
    y_min = -200
    intsteps = 2100

sample, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@SAMPLES_DIR@/dancing.py", script_suffix="@TEST_SUFFIX@",
    cmd_arguments=[sd_method], intsteps=intsteps,
    ref_data="@CMAKE_SOURCE_DIR@/testsuite/python/data/dancing.txt")


@skipIfMissingFeatures
class Sample(ut.TestCase):
    system = sample.system

    def test_trajectory(self):
        # compare the simulated and published trajectories
        simul = sample.positions[:, :, 0:2]
        paper = sample.data.reshape([-1, 3, 2])

        for pid in range(3):
            dist = []
            # the simulated trajectory is oversampled by a ratio of 60:1
            # compared to the published trajectory (60 +/- 5 to 1)
            for desired in paper[:, pid]:
                if desired[1] < y_min:
                    break
                # find the closest point in the simulated trajectory
                idx = np.abs(simul[:, pid, 1] - desired[1]).argmin()
                actual = simul[idx, pid]
                dist.append(np.linalg.norm(actual - desired))
            self.assertLess(idx, sample.intsteps, msg='Insufficient sampling')
            np.testing.assert_allclose(dist, 0, rtol=0, atol=0.7)


if __name__ == "__main__":
    ut.main()
