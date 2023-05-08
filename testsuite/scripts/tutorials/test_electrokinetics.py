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
import importlib_wrapper as iw
import numpy as np

tutorial, skipIfMissingFeatures = iw.configure_and_import(
    "@TUTORIALS_DIR@/electrokinetics/electrokinetics.py", integration_length=400)


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    system = tutorial.system

    def normalize_two_datasets(self, a, b):
        offset = min(np.min(a), np.min(b))
        a -= offset
        b -= offset
        scale = max(np.max(a), np.max(b))
        a /= scale
        b /= scale

    def test_simulation(self):
        for varname, tol in zip(["density", "velocity"], [2, 5]):
            sim = np.array(tutorial.__dict__[varname + "_list"])
            ana = np.array(tutorial.eof_analytical.__dict__[varname + "_list"])
            self.normalize_two_datasets(sim, ana)
            accuracy = np.max(np.abs(sim - ana))
            # expecting at most a few percents deviation
            self.assertLess(accuracy, tol / 100.)


if __name__ == "__main__":
    ut.main()
