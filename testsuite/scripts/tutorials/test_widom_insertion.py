#
# Copyright (C) 2022 The ESPResSo project
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
#

import unittest as ut
import importlib_wrapper 
import numpy as np

tutorial, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "/home/thilo/code/espresso2/espresso/testsuite/scripts/tutorials/local_tutorials/widom_insertion/widom_insertion.py", sample_size=50,
    ci_params={"mesh": (12, 12, 12), "cao": 6, "tune": True}
)


@skipIfMissingFeatures
class Tutorial(ut.TestCase):

    def test(self):
        # simulation data points must be close to the Davies curve
        sim = tutorial.excess_chemical_potentials[:, 0]
        ref = tutorial.davies_equation(tutorial.salt_concentrations)
        self.assertLess(np.std(sim - ref), 0.05)


if __name__ == "__main__":
    ut.main()
