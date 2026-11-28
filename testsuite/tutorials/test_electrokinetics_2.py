#
# Copyright (C) 2019-2026 The ESPResSo project
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
import importlib_wrapper as iw
import numpy as np

tutorial, skipIfMissingFeatures = iw.configure_and_import(
    "@TUTORIALS_DIR@/electrokinetics/electrokinetics_part2_electroosmotic_flow.py")


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    system = tutorial.system

    def test_charge_neutrality(self):
        """
        The FFT-based Poisson solver requires a charge-neutral system.
        """
        self.assertAlmostEqual(tutorial.total_charge, 0., delta=1e-10)

    def test_electroosmotic_flow(self):
        np.testing.assert_allclose(
            np.copy(tutorial.density_eof),
            tutorial.analytic_density_eof,
            rtol=0.005)
        np.testing.assert_allclose(
            np.copy(tutorial.velocity_eof),
            tutorial.analytic_velocity_eof,
            rtol=0.05)
        np.testing.assert_allclose(
            np.copy(tutorial.pressure_tensor_eof),
            tutorial.analytic_pressure_tensor_eof,
            rtol=0.05)

    def test_pressure_driven_flow(self):
        # the ion density profile is unaffected by the body force
        np.testing.assert_allclose(
            np.copy(tutorial.density_pressure),
            tutorial.analytic_density_eof,
            rtol=0.005)
        np.testing.assert_allclose(
            np.copy(tutorial.velocity_pressure),
            tutorial.analytic_velocity_pressure,
            rtol=0.05)
        np.testing.assert_allclose(
            np.copy(tutorial.pressure_tensor_pressure),
            tutorial.analytic_pressure_tensor_pressure,
            rtol=0.05)


if __name__ == "__main__":
    ut.main()
