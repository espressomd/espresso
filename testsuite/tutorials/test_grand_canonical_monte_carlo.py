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
import numpy as np

tutorial, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@TUTORIALS_DIR@/grand_canonical_monte_carlo/grand_canonical_monte_carlo.py",
    p3m_params={"mesh": 10, "cao": 6, "r_cut": 8.22})


@skipIfMissingFeatures
class Tutorial(ut.TestCase):

    def test(self):
        ratios = tutorial.c_monomer.magnitude / \
            (2. * tutorial.salt_concentration_si.magnitude)
        ref_xi = tutorial.analytical_solution(ratios)
        sim_xi_minus = tutorial.partition_coefficients_negatives_array
        sim_xi_plus = tutorial.universal_partion_coefficient_positive
        np.testing.assert_allclose(
            sim_xi_minus, sim_xi_plus, rtol=1e-5, atol=1e-5)
        np.testing.assert_allclose(sim_xi_minus / ref_xi, 2., rtol=0., atol=2.)


if __name__ == "__main__":
    ut.main()
