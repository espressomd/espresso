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
    "@TUTORIALS_DIR@/widom_insertion/widom_insertion.py"
)


@skipIfMissingFeatures
class Tutorial(ut.TestCase):

    def test(self):
        np.testing.assert_allclose(
            tutorial.excess_chemical_potentials[:, 0],
            tutorial.davies_equation(tutorial.salt_concentrations),
            rtol=0.15)


if __name__ == "__main__":
    ut.main()
