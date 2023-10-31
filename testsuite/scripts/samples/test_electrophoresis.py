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

np.random.seed(seed=42)

sample, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@SAMPLES_DIR@/electrophoresis.py", N_SAMPLES=400,
    p3m_params={"prefactor": 1., "mesh": [14, 14, 14], "cao": 1, "tune": False,
                "accuracy": 1e-2, "r_cut": 5.3, "alpha": 0.16})


@skipIfMissingFeatures
class Sample(ut.TestCase):
    system = sample.system

    def test_persistence_length(self):
        # These two values differ due to undersampling, they converge
        # to the same value around N_SAMPLES=1000
        self.assertAlmostEqual(sample.persistence_length, 38.2, delta=1)
        self.assertAlmostEqual(sample.persistence_length_obs, 48.6, delta=1)

    def test_mobility(self):
        self.assertAlmostEqual(sample.mu, 1.05, delta=0.02)

    def test_electrophoresis_gradient(self):
        # the force is applied along the x-axis
        com_vel = np.average(sample.COM_v, axis=0)
        self.assertGreater(abs(com_vel[0]), .9)
        self.assertLess(abs(com_vel[1]), .1)
        self.assertLess(abs(com_vel[2]), .1)


if __name__ == "__main__":
    ut.main()
