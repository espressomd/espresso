
# Copyright (C) 2010-2019 The ESPResSo project
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
import unittest_decorators as utx
import numpy as np

import espressomd
import espressomd.lb
from scipy.optimize import curve_fit

AGRID = .5
N_CELLS = 12
TAU = 0.002
SEED = 1
DENS = 2.4
VISC = 1.8
KT = 0.8


class TestLBPressureACF:
    """Tests that the thermalized LB pressure auto correlation function
    is consistent with the chosen viscosity
    """

    system = espressomd.System(box_l=[AGRID * N_CELLS] * 3)

    system.time_step = TAU
    system.cell_system.skin = 0

    def tearDown(self):
        self.system.actors.clear()
        self.system.thermostat.turn_off()

    def test(self):
        # setup
        system = self.system
        lb = self.lb_class(agrid=AGRID, dens=DENS, visc=VISC,
                           tau=TAU, kT=KT, seed=SEED)
        system.actors.add(lb)
        system.thermostat.set_lb(LB_fluid=lb, seed=SEED + 1)

        # Warmup
        system.integrator.run(500)

        # sampling
        p_global = np.zeros((self.steps, 3, 3))
        p_node = np.zeros((self.steps, 3, 3))

        node = lb[0, 0, 0]

        for i in range(self.steps):
            p_node[i] = node.pressure_tensor
            p_global[i] = lb.pressure_tensor

            system.integrator.run(2)

        # Test that <sigma_[i!=j]> ~=0 and sigma_[ij]=sigma_[ji]
        tol_global = 4 / np.sqrt(self.steps)
        tol_node = tol_global * np.sqrt(N_CELLS**3)

        # check single node
        for i in range(3):
            for j in range(i + 1, 3):
                avg_ij = np.average(p_node[:, i, j])
                avg_ji = np.average(p_node[:, j, i])
                self.assertEqual(avg_ij, avg_ji)

                self.assertLess(avg_ij, tol_node)

        # check system-wide pressure
        for i in range(3):
            for j in range(i + 1, 3):
                avg_ij = np.average(p_global[:, i, j])
                avg_ji = np.average(p_global[:, j, i])
                self.assertEqual(avg_ij, avg_ji)

                self.assertLess(avg_ij, tol_global)

        # Check that stress auto correlatin matches dynamic viscosity
        # eta = V/kT integral(stress acf)
        all_viscs = []
        for i in range(3):
            for j in range(i + 1, 3):

                # Calculate acf
                tmp = np.correlate(
                    p_global[:, i, j], p_global[:, i, j], mode="full")
                acf = tmp[len(tmp) // 2:] / self.steps

                # integrate first part numerically, fit exponential to tail
                t_max_fit = 50 * TAU
                ts = np.arange(0, t_max_fit, 2 * TAU)
                numeric_integral = np.trapz(acf[:len(ts)], dx=2 * TAU)

                # fit tail
                def f(x, a, b): return a * np.exp(-b * x)

                (a, b), _ = curve_fit(f, acf[:len(ts)], ts)
                tail = f(ts[-1], a, b) / b

                integral = numeric_integral + tail

                measured_visc = integral * system.volume() / KT

                self.assertAlmostEqual(
                    measured_visc, VISC * DENS, delta=VISC * DENS * .15)
                all_viscs.append(measured_visc)

        # Check average over xy, xz and yz against tighter limit
        self.assertAlmostEqual(np.average(all_viscs),
                               VISC * DENS, delta=VISC * DENS * .07)


# DISABLE CPU TEST UNTIL #3804 IS SOLVED
class TestLBPressureACFCPU(TestLBPressureACF, ut.TestCase):

    def setUp(self):
        self.lb_class = espressomd.lb.LBFluid
        self.steps = 4000


@utx.skipIfMissingGPU()
class TestLBPressureACFGPU(TestLBPressureACF, ut.TestCase):

    def setUp(self):
        self.lb_class = espressomd.lb.LBFluidGPU
        self.steps = 15000


if __name__ == "__main__":
    ut.main()
