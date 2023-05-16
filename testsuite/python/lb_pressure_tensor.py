#
# Copyright (C) 2010-2022 The ESPResSo project
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
import unittest_decorators as utx
import numpy as np

import espressomd
import espressomd.lb
#import scipy.optimize

N_CELLS = 12


class TestLBPressureTensor:
    """
    Test that the thermalized LB pressure auto correlation function
    is consistent with the chosen viscosity
    """

    params = {'tau': 0.002,
              'agrid': 0.5,
              'density': 2.4,
              'kinematic_viscosity': 1.8,
              'kT': 0.8,
              'seed': 2}
    system = espressomd.System(box_l=[params["agrid"] * N_CELLS] * 3)
    system.time_step = params["tau"]
    system.cell_system.skin = 0

    def tearDown(self):
        self.system.actors.clear()
        self.system.thermostat.turn_off()

    def setUp(self):
        # Setup
        system = self.system
        self.lbf = self.lb_class(**self.params, **self.lb_params)
        system.actors.add(self.lbf)
        system.thermostat.set_lb(LB_fluid=self.lbf, seed=42)

        # Warmup
        system.integrator.run(500)

        # Sampling
        self.p_global = np.zeros((self.steps, 3, 3))
        self.p_node0 = np.zeros((self.steps, 3, 3))
        self.p_node1 = np.zeros((self.steps, 3, 3))

        # Define two sample nodes, at the corner and in the center
        node0 = self.lbf[0, 0, 0]
        node1 = self.lbf[3 * [N_CELLS // 2]]

        for i in range(self.steps):
            self.p_node0[i] = node0.pressure_tensor
            self.p_node1[i] = node1.pressure_tensor
            self.p_global[i] = self.lbf.pressure_tensor

            system.integrator.run(2)

    def assert_allclose_matrix(self, x, y, atol_diag, atol_offdiag):
        """Assert that all elements x_ij, y_ij are close with
        different absolute tolerances for on- and off-diagonal elements.

        """
        assert x.shape == y.shape
        n = min(x.shape)
        mask_offdiag = ~np.identity(n, dtype=bool)

        np.testing.assert_allclose(np.diag(x), np.diag(y), atol=atol_diag)
        np.testing.assert_allclose(
            x[mask_offdiag],
            y[mask_offdiag],
            atol=atol_offdiag)

    def test_averages(self):
        # Sound speed for D3Q19 in LB lattice units
        c_s_lb = np.sqrt(1 / 3)
        # And in MD units
        c_s = c_s_lb * self.lbf.agrid / self.system.time_step

        # Test time average of pressure tensor against expectation ...
        # eq. (19) in ladd01a (https://doi.org/10.1023/A:1010414013942):
        # Pi_eq = rho c_s^2 I + rho u * u = rho c_s^2 I + 2 / V (m u^2 / 2),
        # with 3x3-identity matrix I . Equipartition: m u^2 / 2 = kT /2,
        # Pi_eq = rho c_s^2 I + kT / V
        p_avg_expected = np.diag(
            3 * [self.lbf.density * c_s**2 + self.lbf.kT / self.lbf.agrid**3])

        # ... globally,
        self.assert_allclose_matrix(
            np.mean(self.p_global, axis=0),
            p_avg_expected, atol_diag=c_s_lb**2 * 3, atol_offdiag=c_s_lb**2 / 6)

        # ... for two nodes.
        for time_series in [self.p_node0, self.p_node1]:
            self.assert_allclose_matrix(
                np.mean(time_series, axis=0),
                p_avg_expected, atol_diag=c_s_lb**2 * 200, atol_offdiag=c_s_lb**2 * 9)

        # Test that <sigma_[i!=j]> ~=0 and sigma_[ij]==sigma_[ji] ...
        tol_global = 8 / np.sqrt(self.steps)
        tol_node = tol_global * np.sqrt(N_CELLS**3)

        # ... for the two sampled nodes
        for i in range(3):
            for j in range(i + 1, 3):
                avg_node0_ij = np.average(self.p_node0[:, i, j])
                avg_node0_ji = np.average(self.p_node0[:, j, i])
                avg_node1_ij = np.average(self.p_node1[:, i, j])
                avg_node1_ji = np.average(self.p_node1[:, j, i])

                self.assertEqual(avg_node0_ij, avg_node0_ji)
                self.assertEqual(avg_node1_ij, avg_node1_ji)

                self.assertAlmostEqual(avg_node0_ij, 0., delta=tol_node)
                self.assertAlmostEqual(avg_node1_ij, 0., delta=tol_node)

        # ... for the system-wide pressure tensor
        for i in range(3):
            for j in range(i + 1, 3):
                avg_ij = np.average(self.p_global[:, i, j])
                avg_ji = np.average(self.p_global[:, j, i])
                self.assertEqual(avg_ij, avg_ji)

                self.assertAlmostEqual(avg_ij, 0., delta=tol_node)

        node = self.lbf[0, 0, 0]
        p_eq = np.diag(3 * [self.lbf.density * c_s**2])
        np.testing.assert_allclose(
            np.copy(node.pressure_tensor) - p_eq,
            np.copy(node.pressure_tensor_neq), rtol=0., atol=1e-7)


@utx.skipIfMissingFeatures("WALBERLA")
class TestLBPressureTensorCPU(TestLBPressureTensor, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": True}
    steps = 5000


# TODO WALBERLA
"""
@utx.skipIfMissingFeatures("WALBERLA")
@utx.skipIfMissingGPU()
class TestLBPressureTensorGPU(TestLBPressureTensor, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberlaGPU
    lb_params = {"single_precision": True}
    steps = 50000

    def test_gk_viscosity(self):
        # Check that stress auto correlation matches dynamic viscosity
        # eta = V/kT integral (stress acf), e.g., eq. (5) in Cui et. et al
        # (https://doi.org/10.1080/00268979609484542).
        # Cannot be run for CPU with sufficient statistics without CI timeout.
        dyn_visc = self.params["kinematic_viscosity"] * self.params["density"]
        tau = self.params["tau"]
        kT = self.params["kT"]
        all_viscs = []
        for i in range(3):
            for j in range(i + 1, 3):

                # Calculate acf
                tmp = np.correlate(
                    self.p_global[:, i, j],
                    self.p_global[:, i, j], mode="full")
                acf = tmp[len(tmp) // 2:] / self.steps

                # integrate first part numerically, fit exponential to tail
                t_max_fit = 50 * tau
                ts = np.arange(0, t_max_fit, 2 * tau)
                numeric_integral = np.trapz(acf[:len(ts)], dx=2 * self.params["tau"])

                # fit tail
                def f(x, a, b): return a * np.exp(-b * x)

                (a, b), _ = scipy.optimize.curve_fit(f, acf[:len(ts)], ts)
                tail = f(ts[-1], a, b) / b

                integral = numeric_integral + tail

                measured_visc = integral * self.system.volume() / kT

                self.assertAlmostEqual(
                    measured_visc, dyn_visc, delta=dyn_visc * .15)
                all_viscs.append(measured_visc)

        # Check average over xy, xz and yz against tighter limit
        self.assertAlmostEqual(np.average(all_viscs),
                               dyn_visc, delta=dyn_visc * .07)
"""


if __name__ == "__main__":
    ut.main()
