#
# Copyright (C) 2021-2023 The ESPResSo project
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
import scipy.optimize

sample, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@SAMPLES_DIR@/lb_circular_couette.py")


def taylor_couette(v1, v2, r1, r2, agrid):
    # Taylor-Couette equation
    mu = v2 / v1
    eta = r1 / r2
    scale = 1. / 2**3 / agrid
    a = scale * v1 * (mu - eta**2) / (1 - eta**2)
    b = scale * v1 * r1**2 * (1 - mu) / (1 - eta**2)
    return a, b


@skipIfMissingFeatures
class Sample(ut.TestCase):
    system = sample.system

    def test_taylor_couette_flow(self):
        # get flow profile
        v_r, v_phi, v_z = sample.profile_v.T

        # check velocity is zero for the radial and axial components
        np.testing.assert_allclose(v_r, 0., atol=1e-4)
        np.testing.assert_allclose(v_z, 0., atol=1e-6)

        # check azimuthal velocity is zero inside boundary
        np.testing.assert_allclose(v_phi[:7], 0., atol=1e-7)

        # check azimuthal velocity in the linear regime
        self.assertGreater(v_phi[7], v_phi[6])
        self.assertGreater(v_phi[8], v_phi[7])
        self.assertGreater(v_phi[9], v_phi[8])

        # check azimuthal velocity in the Couette regime
        xdata = sample.profile_r[9:]
        ydata = v_phi[9:]
        a_ref, b_ref = taylor_couette(
            sample.velocity_magnitude, 0.0, sample.cylinder_in.radius,
            sample.cylinder_out.radius, sample.agrid)
        (a_sim, b_sim), _ = scipy.optimize.curve_fit(
            lambda x, a, b: a * x + b / x, xdata, ydata)
        np.testing.assert_allclose([a_sim, b_sim], [a_ref, b_ref], atol=1e-3)


if __name__ == "__main__":
    ut.main()
