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
    "@SAMPLES_DIR@/lb_profile.py", gpu=True)


@skipIfMissingFeatures
class Sample(ut.TestCase):
    system = sample.system

    def test_fit(self):
        # compare the simulated data against the analytical solution
        sim = sample.lb_fluid_profile[:, 0, 0, 2]
        ana = sample.poiseuille_flow(sample.r, sample.r_max, 0.15)
        self.assertLess(np.max(np.abs(sim - ana)), 0.16)


if __name__ == "__main__":
    ut.main()
