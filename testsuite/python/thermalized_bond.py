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


@utx.skipIfMissingFeatures(["MASS"])
class ThermalizedBond(ut.TestCase, ThermostatsCommon):

    """Tests the two velocity distributions for COM and distance created by the
       thermalized bond independently against the single component Maxwell
       distribution. Adapted from langevin_thermostat testcase."""

    box_l = 10.0
    system = espressomd.System(box_l=[box_l] * 3)
    system.cell_system.set_n_square()
    system.cell_system.skin = 0.3

    @classmethod
    def setUpClass(cls):
        np.random.seed(42)

    def test_com_langevin(self):
        """Test for COM thermalization."""

        N = 200
        N2 = int(N / 2)
        self.system.part.clear()
        self.system.time_step = 0.02
        self.system.periodicity = [0, 0, 0]

        m1 = 1.0
        m2 = 10.0
        # Place particles
        for i in range(0, N, 2):
            self.system.part.add(pos=np.random.random(3), mass=m1)
            self.system.part.add(pos=np.random.random(3), mass=m2)

        t_dist = 0
        g_dist = 0
        t_com = 2.0
        g_com = 4.0

        thermalized_dist_bond = espressomd.interactions.ThermalizedBond(
            temp_com=t_com, gamma_com=g_com, temp_distance=t_dist,
            gamma_distance=g_dist, r_cut=2.0, seed=55)
        self.system.bonded_inter.add(thermalized_dist_bond)

        for p1, p2 in zip(self.system.part[::2], self.system.part[1::2]):
            p1.add_bond((thermalized_dist_bond, p2))

        # Warmup
        self.system.integrator.run(50)

        # Sampling
        loops = 150
        v_stored = np.zeros((N2 * loops, 3))
        for i in range(loops):
            self.system.integrator.run(5)
            v_com = 1.0 / \
                (m1 + m2) * \
                (m1 * self.system.part[::2].v + m2 * self.system.part[1::2].v)
            v_stored[i * N2:(i + 1) * N2, :] = v_com

        v_minmax = 5
        bins = 50
        error_tol = 0.017
        self.check_velocity_distribution(
            v_stored, v_minmax, bins, error_tol, t_com)

    def test_dist_langevin(self):
        """Test for dist thermalization."""

        N = 100
        N2 = int(N / 2)
        self.system.part.clear()
        self.system.time_step = 0.02
        self.system.periodicity = [1, 1, 1]

        m1 = 1.0
        m2 = 10.0
        # Place particles
        for i in range(0, N, 2):
            self.system.part.add(pos=np.random.random(3), mass=m1)
            self.system.part.add(pos=np.random.random(3), mass=m2)

        t_dist = 2.0
        g_dist = 4.0
        t_com = 0.0
        g_com = 0.0

        thermalized_dist_bond = espressomd.interactions.ThermalizedBond(
            temp_com=t_com, gamma_com=g_com, temp_distance=t_dist,
            gamma_distance=g_dist, r_cut=9, seed=51)
        self.system.bonded_inter.add(thermalized_dist_bond)

        for p1, p2 in zip(self.system.part[::2], self.system.part[1::2]):
            p1.add_bond((thermalized_dist_bond, p2))

        # Warmup
        self.system.integrator.run(50)

        # Sampling
        loops = 150
        v_stored = np.zeros((N2 * loops, 3))
        for i in range(loops):
            self.system.integrator.run(5)
            v_dist = self.system.part[1::2].v - self.system.part[::2].v
            v_stored[i * N2:(i + 1) * N2, :] = v_dist

        v_minmax = 5
        bins = 50
        error_tol = 0.015
        self.check_velocity_distribution(
            v_stored, v_minmax, bins, error_tol, t_dist)


if __name__ == "__main__":
    ut.main()
