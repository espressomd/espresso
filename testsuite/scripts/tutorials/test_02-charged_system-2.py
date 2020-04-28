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
    "@TUTORIALS_DIR@/02-charged_system/02-charged_system-2.py",
    num_steps_equilibration=200, num_configs=5, integ_steps_per_config=60)


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    system = tutorial.system

    def test_distribution(self):
        """
        checks if the particle distribution is within the box
        """
        for i in range(1, 3):
            pos = np.flatnonzero(tutorial.res[:, i] > 0)
            self.assertGreater(tutorial.res[pos[0], 0], tutorial.wall_margin)
            self.assertLess(tutorial.res[pos[-1], 0],
                            tutorial.box_z - tutorial.wall_margin)


if __name__ == "__main__":
    ut.main()
