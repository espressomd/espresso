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
import scipy.optimize

tutorial, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@TUTORIALS_DIR@/lattice_boltzmann/lattice_boltzmann_part2.py")


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    system = tutorial.system

    def test_ballistic_regime(self):
        for (tau_p, tau, msd) in zip(tutorial.tau_p_values,
                                     tutorial.tau_results,
                                     tutorial.msd_results):
            popt, _ = scipy.optimize.curve_fit(
                tutorial.quadratic, tau[:tau_p], msd[:tau_p])
            residuals = msd[:tau_p] - tutorial.quadratic(tau[:tau_p], *popt)
            np.testing.assert_allclose(residuals, 0, rtol=0, atol=1e-3)

    def test_diffusion_coefficient(self):
        D_val = tutorial.diffusion_results
        D_ref = tutorial.KT / np.array(tutorial.gammas)
        np.testing.assert_allclose(D_val, D_ref, rtol=0, atol=0.1)


if __name__ == "__main__":
    ut.main()
