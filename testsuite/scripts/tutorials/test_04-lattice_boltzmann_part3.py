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

tutorial, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@TUTORIALS_DIR@/04-lattice_boltzmann/04-lattice_boltzmann_part3.py",
    LOOPS=2000, gpu=True)


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


if __name__ == "__main__":
    ut.main()
