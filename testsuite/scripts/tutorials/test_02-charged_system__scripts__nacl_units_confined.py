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
    "@TUTORIALS_DIR@/02-charged_system/scripts/nacl_units_confined.py",
    num_steps_equilibration=1000, num_configs=40, integ_steps_per_config=60,
    l_bjerrum=130878 / 2.)


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    system = tutorial.system

    def test_distance_Na_Cl(self):
        d_Na_Electrode = np.argmax(tutorial.z_dens_na)
        d_Cl_Electrode = np.argmax(tutorial.z_dens_cl)
        self.assertGreaterEqual(d_Na_Electrode, 88)
        self.assertLessEqual(d_Cl_Electrode, 12)


if __name__ == "__main__":
    ut.main()
