#
# Copyright (C) 2022 The ESPResSo project
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
import importlib_wrapper
import numpy as np
import scipy.stats


tutorial, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@TUTORIALS_DIR@/lattice_boltzmann/lattice_boltzmann_sedimentation.py",
    sampling_steps=450)


def curl2d(flow, spacing):
    # curl_z = d(v_y) / dx - d(v_x) / dy
    dvy_dx = np.gradient(flow[:, :, 1], spacing[0], axis=0)
    dvx_dy = np.gradient(flow[:, :, 0], spacing[1], axis=1)
    return dvy_dx - dvx_dy


def get_peak_position(mat, kernel):
    return np.unravel_index(kernel(mat, axis=None), mat.shape)


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    system = tutorial.system

    def test_flow_profile(self):
        # slice trajectory to keep only the flow field onset
        flow_field = tutorial.data_flowfield[400:450, :, :, 0:2]
        vortices = np.zeros((2, flow_field.shape[0], 2))
        for i in range(flow_field.shape[0]):
            curl = curl2d(flow_field[i], 2 * [tutorial.spacing])
            vortices[0, i] = get_peak_position(curl, np.argmax)
            vortices[1, i] = get_peak_position(curl, np.argmin)

        # check flow field curl
        ref_pos_x = [5.5, 13.5]  # LB units
        ref_pos_y = 2 * [13.5]   # LB units
        for i in range(2):
            width = tutorial.n_width  # LB units
            vortex_avg_x = scipy.stats.circmean(vortices[i, :, 0], high=width)
            vortex_std_x = scipy.stats.circstd(vortices[i, :, 0], high=width)
            vortex_avg_y = np.mean(vortices[i, :, 1])
            vortex_std_y = np.std(vortices[i, :, 1])
            self.assertAlmostEqual(vortex_avg_x, ref_pos_x[i], delta=2.)
            self.assertAlmostEqual(vortex_avg_y, ref_pos_y[i], delta=2.)
            self.assertLess(vortex_std_x, 2.)
            self.assertLess(vortex_std_y, 3.)


if __name__ == "__main__":
    ut.main()
