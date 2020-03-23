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
    "@TUTORIALS_DIR@/01-lennard_jones/01-lennard_jones.py",
    substitutions=lambda code: code.replace('r_star = ', 'r_star = r #', 1))


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    system = tutorial.system

    def test_global_variables(self):
        self.assertLess(tutorial.standard_error_total_energy, 2.5)

    def test_rdf_curve(self):
        self.assertGreater(
            np.corrcoef(tutorial.avg_rdf, tutorial.theo_rdf)[1, 0], 0.985)


if __name__ == "__main__":
    ut.main()
