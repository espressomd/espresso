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

tutorial, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@TUTORIALS_DIR@/02-charged_system/02-charged_system-1.py",
    N_frames=20, steps_per_frame=100, warmup_steps=1000,
    N_frames_salt=200)


@skipIfMissingFeatures
class Tutorial02_1(ut.TestCase):
    system = tutorial.system

    def test_overcharging(self):
        """ 
        Test that adding salt leads to a positive layer around the rod
        """

        self.assertGreater(max(tutorial.charge_hist), 1.05)


if __name__ == "__main__":
    ut.main()
