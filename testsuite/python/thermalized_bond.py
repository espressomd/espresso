#
# Copyright (C) 2013-2022 The ESPResSo project
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
import thermostats_common

import espressomd
import espressomd.interactions


@utx.skipIfMissingFeatures(["MASS"])
class ThermalizedBond(ut.TestCase, thermostats_common.ThermostatsCommon):

    """
    Test the two velocity distributions for COM and distance created by the
    thermalized bond independently against the single component Maxwell
    distribution.
    """

    system = espressomd.System(box_l=[10., 10., 10.])
    system.cell_system.set_n_square()
    system.cell_system.skin = 0.3

    def setUp(self):
        np.random.seed(42)
        self.system.time_step = 0.01
        self.system.periodicity = [True, True, True]

    def tearDown(self):
        self.system.part.clear()
        self.system.thermostat.turn_off()

    def test_com_langevin(self):
        """Test for COM thermalization."""

        N = 200
        N_half = int(N / 2)
        self.system.time_step = 0.02
        self.system.periodicity = [False, False, False]

        m1 = 1.0
        m2 = 10.0
        # Place particles

        partcls_m1 = self.system.part.add(
            pos=np.random.random(
                (N_half, 3)), mass=N_half * [m1])
        partcls_m2 = self.system.part.add(
            pos=np.random.random(
                (N_half, 3)), mass=N_half * [m2])

        t_dist = 0
        g_dist = 0
        t_com = 2.0
        g_com = 4.0

        self.system.thermostat.set_thermalized_bond(seed=55)
        thermalized_dist_bond = espressomd.interactions.ThermalizedBond(
            temp_com=t_com, gamma_com=g_com, temp_distance=t_dist,
            gamma_distance=g_dist, r_cut=2.0)
        self.system.bonded_inter.add(thermalized_dist_bond)

        for p1, p2 in zip(partcls_m1, partcls_m2):
            p1.add_bond((thermalized_dist_bond, p2))

        # Warmup
        self.system.integrator.run(50)

        # Sampling
        loops = 150
        v_stored = np.zeros((N_half * loops, 3))
        for i in range(loops):
            self.system.integrator.run(5)
            v_com = 1.0 / \
                (m1 + m2) * \
                (m1 * partcls_m1.v + m2 * partcls_m2.v)
            v_stored[i * N_half:(i + 1) * N_half, :] = v_com

        v_minmax = 5
        bins = 50
        error_tol = 0.017
        self.check_velocity_distribution(
            v_stored, v_minmax, bins, error_tol, t_com)

    def test_dist_langevin(self):
        """Test for dist thermalization."""

        N = 100
        N_half = int(N / 2)
        self.system.time_step = 0.02
        self.system.periodicity = [True, True, True]

        m1 = 1.0
        m2 = 10.0
        # Place particles
        partcls_m1 = self.system.part.add(
            pos=np.random.random(
                (N_half, 3)), mass=N_half * [m1])
        partcls_m2 = self.system.part.add(
            pos=np.random.random(
                (N_half, 3)), mass=N_half * [m2])

        t_dist = 2.0
        g_dist = 4.0
        t_com = 0.0
        g_com = 0.0

        self.system.thermostat.set_thermalized_bond(seed=51)
        thermalized_dist_bond = espressomd.interactions.ThermalizedBond(
            temp_com=t_com, gamma_com=g_com, temp_distance=t_dist,
            gamma_distance=g_dist, r_cut=9)
        self.system.bonded_inter.add(thermalized_dist_bond)

        for p1, p2 in zip(partcls_m1, partcls_m2):
            p1.add_bond((thermalized_dist_bond, p2))

        # Warmup
        self.system.integrator.run(50)

        # Sampling
        loops = 150
        v_stored = np.zeros((N_half * loops, 3))
        for i in range(loops):
            self.system.integrator.run(5)
            v_dist = partcls_m2.v - partcls_m1.v
            v_stored[i * N_half:(i + 1) * N_half, :] = v_dist

        v_minmax = 5
        bins = 50
        error_tol = 0.015
        self.check_velocity_distribution(
            v_stored, v_minmax, bins, error_tol, t_dist)

    def test_exceptions(self):
        espressomd.interactions.ThermalizedBond(
            temp_com=1., gamma_com=1., temp_distance=1., gamma_distance=1.,
            r_cut=3.)
        with self.assertRaisesRegex(Exception, "Thermalized bonds require the thermalized_bond thermostat"):
            self.system.integrator.run(0)


if __name__ == "__main__":
    ut.main()
