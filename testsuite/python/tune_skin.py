# Copyright (C) 2010-2019 The ESPResSo project
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
        with self.assertRaisesRegex(Exception, 'Error during tune_skin'):
            self.system.cell_system.tune_skin(
                min_skin=0.1,
                max_skin=0.6,
                tol=0.05,
                int_steps=3)

    def test_works_with_adjustment(self):
        self.system.cell_system.tune_skin(
            min_skin=0.1,
            max_skin=0.6,
            tol=0.05,
            int_steps=3,
            adjust_max_skin=True)


if __name__ == "__main__":
    ut.main()
