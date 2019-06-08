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
import importlib_wrapper as iw
import numpy as np
import sys

# these tutorials need to be executed sequentially
tutorial_simulation, skipIfMissingFeatures_simulation = iw.configure_and_import(
    "@TUTORIALS_DIR@/07-electrokinetics/scripts/eof_electrokinetics.py",
    gpu=True, integration_length=600, dt=0.5)
# use importlib directly to avoid an error for some myconfig.hpp configurations
sys.path.insert(0, "@TUTORIALS_DIR@/07-electrokinetics/scripts/")
tutorial_analytical = iw.importlib.import_module("eof_analytical")
tutorial_plot, skipIfMissingFeatures_plot = iw.configure_and_import(
    "@TUTORIALS_DIR@/07-electrokinetics/scripts/plot.py")


@skipIfMissingFeatures_simulation
@skipIfMissingFeatures_plot
class Tutorial(ut.TestCase):
    system = tutorial_simulation.system

    def normalize_two_datasets(self, a, b):
        offset = min(np.min(a), np.min(b))
        a -= offset
        b -= offset
        scale = max(np.max(a), np.max(b))
        a /= scale
        b /= scale

    def test_simulation(self):
        for varname in ("density", "velocity", "pressure_xz"):
            sim = np.array(tutorial_simulation.__dict__[varname + "_list"])
            ana = np.array(tutorial_analytical.__dict__[varname + "_list"])
            self.normalize_two_datasets(sim, ana)
            accuracy = np.max(np.abs(sim - ana))
            # expecting at most 3% deviation
            self.assertLess(accuracy, 3.0 / 100)


if __name__ == "__main__":
    ut.main()
