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
    "@TUTORIALS_DIR@/11-ferrofluid/11-ferrofluid_part3.py",
    equil_steps=200,equil_rounds=10,loops=250,alphas=[0.5])


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    system = tutorial.system


    def test(self):
        self.assertGreater(
            tutorial.magnetization_star[1],
            tutorial.L(tutorial.alpha_mean_field(tutorial.alphas[1],tutorial.dip_lambda,tutorial.phi)))
        self.assertLess(
            tutorial.magnetization_star[1],
            1)
        self.assertAlmostEqual(
            tutorial.chi,tutorial.chi_I,delta=0.05)
        self.assertGreater(
            tutorial.chi,tutorial.chi_L)



if __name__ == "__main__":
    ut.main()
