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

if '@TEST_SUFFIX@' == 'rouse':
    params = {}
elif '@TEST_SUFFIX@' == 'zimm':
    params = {'LOOPS': 400, 'POLYMER_MODEL': 'Zimm'}

tutorial, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@TUTORIALS_DIR@/polymers/polymers.py",
    script_suffix="@TEST_SUFFIX@", **params)


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    system = tutorial.system

    def test_exponents(self):
        msg = 'The R_F exponent should be close to 0.588'
        self.assertGreater(tutorial.rf_exponent, 0.50, msg=msg)
        self.assertLess(tutorial.rf_exponent, 0.85, msg=msg)
        msg = 'The R_g exponent should be close to 0.588'
        self.assertGreater(tutorial.rg_exponent, 0.50, msg=msg)
        self.assertLess(tutorial.rg_exponent, 0.75, msg=msg)
        msg = 'The R_h exponent should be close to 0.333'
        self.assertGreater(tutorial.rh_exponent, 0.30, msg=msg)
        self.assertLess(tutorial.rh_exponent, 0.50, msg=msg)
        np.testing.assert_allclose(tutorial.rf2_rg2_ratio, 6.0, atol=1.1,
                                   err_msg='R_F^2/R_g^2 should be close to 6.0')

    def test_diffusion_coefficients(self):
        # polymer diffusion
        ref_D = [0.0363, 0.0269, 0.0234]
        np.testing.assert_allclose(tutorial.diffusion_msd, ref_D, rtol=0.30)
        np.testing.assert_allclose(tutorial.diffusion_gk, ref_D, rtol=0.15)
        # monomer diffusion
        if tutorial.POLYMER_MODEL == 'Rouse':
            ref_D0 = tutorial.KT / tutorial.GAMMA
            self.assertAlmostEqual(tutorial.popt_msd[0], ref_D0, delta=0.02)
            self.assertAlmostEqual(tutorial.popt_gk[0], ref_D0, delta=0.02)


if __name__ == "__main__":
    ut.main()
