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

sample, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@SAMPLES_DIR@/electrophoresis.py")


@skipIfMissingFeatures
class Sample(ut.TestCase):
    system = sample.system

    @ut.skipIf(tuple(map(int, np.__version__.split("."))) < (1, 10),
               "not supported for numpy < v1.10")
    def test_persistence_length(self):
        value = np.mean(sample.fit[0])
        self.assertTrue(5. < value < 110., "length = {:.0f}".format(value))

    def test_mobility(self):
        value = sample.mu
        self.assertTrue(-0.4 < value < 0.4, "mobility = {:.2f}".format(value))

    def test_electrophoresis_gradient(self):
        # the force is applied along the x-axis
        gradient = np.mean(np.gradient(sample.COM.T, axis=1), axis=1)
        self.assertLess(abs(gradient[0] + 1e-2), 2e-3)
        self.assertLess(abs(gradient[1]), 2e-3)
        self.assertLess(abs(gradient[2]), 2e-3)


if __name__ == "__main__":
    ut.main()
