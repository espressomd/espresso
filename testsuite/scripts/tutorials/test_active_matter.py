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

import numpy as np
import unittest as ut
import importlib_wrapper
import os

tutorial, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@TUTORIALS_DIR@/active_matter/active_matter.py",
    ED_N_SAMPLING_STEPS=100000,
    RECT_N_SAMPLES=150,
    HYDRO_N_STEPS=150
)


def curl2d(flow, spacing):
    # curl_z = d(v_y) / dx - d(v_x) / dy
    dvy_dx = np.gradient(flow[:, :, 1], spacing[0], axis=0)
    dvx_dy = np.gradient(flow[:, :, 0], spacing[1], axis=1)
    return dvy_dx - dvx_dy


@skipIfMissingFeatures
class TestActMat(ut.TestCase):
    system = tutorial.system

    def test_enhanced_diffusion(self):
        """ Check that the active particle diffuses faster than the passive one
        """
        self.assertGreater(
            tutorial.msd_result[-1, 0], tutorial.msd_result[-1, 1])

    def test_rectification(self):
        """ Check that the center of mass is in the right half of the box
        """
        self.assertGreater(tutorial.com_deviations[-1], 0)

    def test_hydrodynamics(self):
        """ Check that the particle is moving up and the fluid down
        """
        self.assertGreater(
            tutorial.system.analysis.linear_momentum(
                include_lbfluid=False)[2], 0)
        self.assertLess(
            tutorial.system.analysis.linear_momentum(
                include_particles=False)[2], 0)

    def test_flow_profile(self):
        """
        Check the flow field curl is symmetric with 4 local extrema
        centered around the particle.
        """
        flow_field = tutorial.vels[:, :, (0, 2)]
        curl = curl2d(flow_field, 2 * [tutorial.lbf.agrid])
        curl_percent = curl * 100. / np.max(np.abs(curl))
        threshold_percent = 85.
        self.assertGreaterEqual(curl_percent[16, 20], threshold_percent)
        self.assertGreaterEqual(curl_percent[18, 16], threshold_percent)
        self.assertLessEqual(curl_percent[16, 16], -threshold_percent)
        self.assertLessEqual(curl_percent[18, 20], -threshold_percent)

    def test_file_generation(self):
        for name in ["position_0.vtk", "lb_velocity_0.vtu"]:
            filepath = os.path.join(tutorial.vtk_outdir, name)
            self.assertTrue(
                os.path.isfile(filepath),
                filepath + " not created")


if __name__ == "__main__":
    ut.main()
