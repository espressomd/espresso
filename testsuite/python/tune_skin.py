#
# Copyright (C) 2010-2026 The ESPResSo project
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
import unittest_decorators as utx
import espressomd


@utx.skipIfMissingFeatures("LENNARD_JONES")
class TuneSkin(ut.TestCase):
    system = espressomd.System(box_l=[1.35, 2.4, 1.7])
    system.time_step = 0.01

    def setUp(self):
        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=1,
            sigma=0.2,
            cutoff=0.3,
            shift="auto")

    def test_fails_without_adjustment(self):
        with self.assertRaisesRegex(Exception, "number of cells 2 is smaller than minimum 8: either interaction range is too large for the current skin .+ or min_num_cells too large"):
            self.system.cell_system.tune_skin(
                min_skin=0.1,
                max_skin=0.6,
                tol=0.05,
                int_steps=3)

    def test_works_with_adjustment(self):
        skin = self.system.cell_system.tune_skin(
            min_skin=0.1,
            max_skin=0.6,
            tol=0.05,
            int_steps=3,
            adjust_max_skin=True)
        self.assertAlmostEqual(skin, self.system.cell_system.skin, delta=1e-12)

    def test_rejects_invalid_parameters(self):
        # A non-positive tolerance must be rejected with an exception. Without
        # validation the core bisection loop `while (fabs(a - b) > tol)` never
        # terminates for tol <= 0 (fabs(a - b) >= 0 is always > a negative tol,
        # and once the interval collapses the midpoint stops progressing), so
        # the call would hang forever instead of returning.
        for tol in (-1.0, 0.0):
            with self.assertRaisesRegex(ValueError, "Parameter 'tol' must be > 0"):
                self.system.cell_system.tune_skin(
                    min_skin=0.1,
                    max_skin=0.6,
                    tol=tol,
                    int_steps=3,
                    adjust_max_skin=True)
        # A non-positive integration-step count is meaningless and leads to a
        # division by zero in the per-step timing, so it must be rejected.
        for int_steps in (-1, 0):
            with self.assertRaisesRegex(ValueError, "Parameter 'int_steps' must be > 0"):
                self.system.cell_system.tune_skin(
                    min_skin=0.1,
                    max_skin=0.6,
                    tol=0.05,
                    int_steps=int_steps,
                    adjust_max_skin=True)
        with self.assertRaisesRegex(ValueError, "Parameter 'min_skin' must be >= 0"):
            self.system.cell_system.tune_skin(
                min_skin=-0.1, max_skin=0.6, tol=0.01, int_steps=3)
        with self.assertRaisesRegex(ValueError, "Parameter 'max_skin' must be >= 0"):
            self.system.cell_system.tune_skin(
                min_skin=0.1, max_skin=-0.6, tol=0.01, int_steps=3)
        with self.assertRaisesRegex(ValueError, "Parameter 'max_skin' must be >= 'min_skin'"):
            self.system.cell_system.tune_skin(
                min_skin=0.6, max_skin=0.1, tol=0.01, int_steps=3)


if __name__ == "__main__":
    ut.main()
