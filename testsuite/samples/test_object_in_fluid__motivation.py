# Copyright (C) 2019-2022 The ESPResSo project
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
import pathlib
import numpy as np

os.chdir("@SAMPLES_DIR@/object_in_fluid")
sample, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@SAMPLES_DIR@/object_in_fluid/motivation.py", maxCycle=4, cmd_arguments=[0])


@skipIfMissingFeatures
class Sample(ut.TestCase):
    system = sample.system

    def test_file_generation(self):
        basenames = [
            "cylinderA.vtk",
            "cylinderB.vtk",
            "cylinderC.vtk",
            "wallTop.vtk",
            "wallBottom.vtk",
            "wallBack.vtk",
            "wallFront.vtk"]
        for j in range(sample.maxCycle + 1):
            basenames.append(f"cell0_{j}.vtk")

        # test .vtk files exist
        path_vtk_root = pathlib.Path("output")
        for name in basenames:
            filepath = path_vtk_root / "sim0" / name
            self.assertTrue(
                filepath.is_file(),
                f"File '{filepath}' not created")

        # make sure we are still in the LB linear regime
        lb_density = np.copy(sample.lbf[:, :, :].density)
        np.testing.assert_allclose(lb_density, 1., atol=2e-3)

        # verify cell momentum
        cell_vel = np.mean(self.system.part.all().v, axis=0)
        np.testing.assert_allclose(
            cell_vel, [1.54e-2, 1.8e-3, 0.], rtol=1e-2, atol=1e-6)


if __name__ == "__main__":
    ut.main()
