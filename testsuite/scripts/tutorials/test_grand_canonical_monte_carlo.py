#
# Copyright (C) 2023 The ESPResSo project
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
import pint
import numpy as np

ureg = pint.UnitRegistry()
tutorial, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@TUTORIALS_DIR@/grand_canonical_monte_carlo/grand_canonical_monte_carlo.py",
    salt_concentration_magnitudes_si=[0.003, 0.01], number_of_loops_for_each_concentration=[250, 100],
    excess_chemical_potential_data=[-0.123732052028611, -0.218687259792629],
    excess_chemical_potential_data_error=[0.00152160176511698, 0.00220162667953136])


@skipIfMissingFeatures
class Tutorial(ut.TestCase):

    def test(self):
        sim_xi_minus = tutorial.partition_coefficients_negatives_array
        self.assertLess(np.abs(sim_xi_minus[0] - 0.05), 0.02)
        self.assertLess(np.abs(sim_xi_minus[1] - 0.15), 0.1)


if __name__ == "__main__":
    ut.main()
