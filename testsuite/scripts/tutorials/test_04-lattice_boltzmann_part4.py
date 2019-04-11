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


def substitutions(code):
    # inject the solution into the notebook
    with open("@TUTORIALS_DIR@/04-lattice_boltzmann/scripts/part4_plot.py") as f:
        code_plot = f.read()
    return code + "\n" + code_plot


tutorial, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@TUTORIALS_DIR@/04-lattice_boltzmann/04-lattice_boltzmann_part4.py",
    gpu=True, substitutions=substitutions)


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    system = tutorial.system

    def test_analytical_solution(self):
        difference = tutorial.scaling * tutorial.sim - tutorial.ana
        accuracy = np.max(np.abs(difference))
        self.assertLess(accuracy, 1e-3)


if __name__ == "__main__":
    ut.main()
