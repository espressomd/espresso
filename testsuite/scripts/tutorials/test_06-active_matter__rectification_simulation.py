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
import os
import numpy as np

tutorial, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@TUTORIALS_DIR@/06-active_matter/SOLUTIONS/rectification_simulation.py",
    cmd_arguments=[6.0], prod_steps=80, prod_length=100)


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    system = tutorial.system

    def test_rectification(self):
        x = tutorial.system.part[:].pos[:, 0]
        left_chamber = np.sum(x < tutorial.length / 2.0)
        right_chamber = np.sum(x > tutorial.length / 2.0)
        excess = (right_chamber - left_chamber) * 100. / tutorial.npart
        # expecting at least 5% excess due to rectification
        self.assertGreater(excess, 5.0)

    def test_file_generation(self):
        # test .vtk/.dat files exist
        for name in ["CMS_{}.dat", "points_{}.vtk"]:
            filepath = os.path.join(tutorial.outdir, name.format(tutorial.vel))
            self.assertTrue(
                os.path.isfile(filepath),
                filepath + " not created")


if __name__ == "__main__":
    ut.main()
