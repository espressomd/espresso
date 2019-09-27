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
    "@TUTORIALS_DIR@/04-lattice_boltzmann/scripts/04-lattice_boltzmann_part4_solution.py",
    gpu=True)


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    def test_flow_profile(self):
        analytical = tutorial.poiseuille_flow(tutorial.x_values - (tutorial.HEIGHT / 2.0 + tutorial.AGRID),
                                              tutorial.FORCE_DENSITY[1],
                                              tutorial.VISCOSITY * tutorial.DENSITY,
                                              tutorial.HEIGHT)
        simulation = tutorial.fluid_velocities[:, 1]
        rmsd = np.sqrt(np.sum(np.square(analytical - simulation)))
        self.assertLess(rmsd, 0.02 * tutorial.AGRID / tutorial.lbf.tau)


if __name__ == "__main__":
    ut.main()
