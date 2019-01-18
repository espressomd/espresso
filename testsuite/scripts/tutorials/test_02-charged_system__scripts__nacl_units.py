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
    "@TUTORIALS_DIR@/02-charged_system/scripts/nacl_units.py",
    num_steps_equilibration=100, num_configs=20, integ_steps_per_config=80)


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    system = tutorial.system

    def test_crystal(self):
        lattice_constant = 2.63
        d_Na_Cl = tutorial.r[np.argmax(tutorial.rdf_01)]
        d_Na_Na = tutorial.r[np.argmax(tutorial.rdf_11)]
        d_Cl_Cl = tutorial.r[np.argmax(tutorial.rdf_00)]
        self.assertLess(abs(d_Na_Cl - lattice_constant), 0.1)
        self.assertLess(abs(d_Na_Na - d_Cl_Cl), 0.1)
        self.assertLess(abs(d_Na_Cl - d_Cl_Cl / 2**0.5), 0.1)


if __name__ == "__main__":
    ut.main()
