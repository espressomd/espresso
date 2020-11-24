#
# Copyright (C) 2013-2019 The ESPResSo project
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
import numpy as np
import unittest as ut
import unittest_decorators as utx

import espressomd
from thermostats_common import ThermostatsCommon


@utx.skipIfMissingFeatures("DPD")
class DPDThermostat(ut.TestCase, ThermostatsCommon):

    """Tests the velocity distribution created by the DPD thermostat."""

    s = espressomd.System(box_l=3 * [10.0])
    s.time_step = 0.01
    s.cell_system.skin = 0.4

    def setUp(self):
        np.random.seed(16)

    def tearDown(self):
        s = self.s
        s.part.clear()
        s.thermostat.turn_off()
        s.integrator.set_vv()

    def check_total_zero(self):
        v_total = np.sum(self.s.part[:].v, axis=0)
        np.testing.assert_allclose(v_total, np.zeros(3), atol=1e-11)

    def single(self, with_langevin=False):
        """Test velocity distribution of a dpd fluid with a single type."""
        N = 500
        s = self.s
        s.part.add(pos=s.box_l * np.random.random((N, 3)))
        kT = 2.3
        gamma = 1.5
        if with_langevin:
            s.thermostat.set_langevin(kT=kT, gamma=gamma, seed=41)
        s.thermostat.set_dpd(kT=kT, seed=42)
        s.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=0, gamma=gamma, r_cut=1.5,
            trans_weight_function=0, trans_gamma=gamma, trans_r_cut=1.5)
        s.integrator.run(100)
        loops = 100
        v_stored = np.zeros((loops, N, 3))
        for i in range(loops):
            s.integrator.run(10)
            v_stored[i] = s.part[:].v
        v_minmax = 5
        bins = 5
        error_tol = 0.01
        self.check_velocity_distribution(
            v_stored.reshape((-1, 3)), v_minmax, bins, error_tol, kT)

        if not with_langevin:
            self.check_total_zero()

    def test_single(self):
        self.single()

    def test_single_with_langevin(self):
        self.single(True)

    def test_binary(self):
        """Test velocity distribution of binary dpd fluid"""
        N = 200
        s = self.s
        s.part.add(pos=s.box_l * np.random.random((N // 2, 3)),
                   type=N // 2 * [0])
        s.part.add(pos=s.box_l * np.random.random((N // 2, 3)),
                   type=N // 2 * [1])
        kT = 2.3
        gamma = 3.5
        s.thermostat.set_dpd(kT=kT, seed=42)
        s.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=0, gamma=gamma, r_cut=1.0,
            trans_weight_function=0, trans_gamma=gamma, trans_r_cut=1.0)
        s.non_bonded_inter[1, 1].dpd.set_params(
            weight_function=0, gamma=gamma, r_cut=1.0,
            trans_weight_function=0, trans_gamma=gamma, trans_r_cut=1.0)
        s.non_bonded_inter[0, 1].dpd.set_params(
            weight_function=0, gamma=gamma, r_cut=1.5,
            trans_weight_function=0, trans_gamma=gamma, trans_r_cut=1.5)
        s.integrator.run(100)
        loops = 250
        v_stored = np.zeros((loops, N, 3))
        for i in range(loops):
            s.integrator.run(10)
            v_stored[i] = s.part[:].v
        v_minmax = 5
        bins = 5
        error_tol = 0.01
        self.check_velocity_distribution(
            v_stored.reshape((-1, 3)), v_minmax, bins, error_tol, kT)
        self.check_total_zero()

    def test_disable(self):
        N = 200
        s = self.s
        s.time_step = 0.01
        s.part.add(pos=s.box_l * np.random.random((N, 3)))
        kT = 2.3
        gamma = 1.5
        s.thermostat.set_dpd(kT=kT, seed=42)
        s.non_bonded_inter[0, 0].dpd.set_params(
            weight_function=0, gamma=gamma, r_cut=1.5,
            trans_weight_function=0, trans_gamma=gamma, trans_r_cut=1.5)

        s.integrator.run(10)

        s.thermostat.turn_off()

        # Reset velocities
        s.part[:].v = [1., 2., 3.]

        s.integrator.run(10)

        # Check that there was neither noise nor friction
        for v in s.part[:].v:
            for i in range(3):
                self.assertEqual(v[i], float(i + 1))

        # Turn back on
        s.thermostat.set_dpd(kT=kT, seed=42)

        # Reset velocities for faster convergence
        s.part[:].v = [0., 0., 0.]

        # Equilibrate
        s.integrator.run(250)

        loops = 250
        v_stored = np.zeros((loops, N, 3))
        for i in range(loops):
            s.integrator.run(10)
            v_stored[i] = s.part[:].v
        v_minmax = 5
        bins = 5
        error_tol = 0.012
        self.check_velocity_distribution(
            v_stored.reshape((-1, 3)), v_minmax, bins, error_tol, kT)


if __name__ == "__main__":
    ut.main()
