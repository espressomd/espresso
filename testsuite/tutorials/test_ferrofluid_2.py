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


tutorial, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@TUTORIALS_DIR@/ferrofluid/ferrofluid_part2.py",
    equil_steps=200, equil_rounds=10, loops=500, alphas=[0.5])


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    system = tutorial.system

    def test(self):
        langevin_magnetization_curve = tutorial.L(np.array(tutorial.alphas))
        self.assertGreater(
            tutorial.magnetization_para[0],
            tutorial.magnetization_perp[0])
        self.assertGreater(
            tutorial.magnetization_para[0] / tutorial.N_PART,
            langevin_magnetization_curve)
        self.assertLess(
            tutorial.magnetization_perp[0] / tutorial.N_PART,
            langevin_magnetization_curve)


if __name__ == "__main__":
    ut.main()
