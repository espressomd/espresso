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

import unittest as ut

import numpy as np
import pickle

import espressomd
import espressomd.observables
import espressomd.accumulators


class CorrelatorTest(ut.TestCase):

    """
    Test class for the Correlator accumulator.

    """
    # Handle for espresso system
    system = espressomd.System(box_l=[10, 10, 10])
    system.cell_system.skin = 0.4
    system.time_step = 0.01

    def tearDown(self):
        self.system.part.clear()
        self.system.auto_update_accumulators.clear()

    def calc_tau(self, time_step, tau_lin, length, delta_N=1):
        tau = []
        for i in range(tau_lin):
            tau.append(i)
        factor = 1
        while len(tau) < length:
            p = tau[-1] + factor * 1
            for i in range(0, tau_lin, 2):
                tau.append(p + factor * i)
            factor *= 2
        return time_step * np.array(tau[:length]) * delta_N

    def check_sizes(self, acc, steps, linear=False):
        sizes = acc.sample_sizes()
        tau_lin = acc.tau_lin
        max_lin = np.arange(steps, steps - tau_lin - 1, -1)
        if linear:
            np.testing.assert_equal(sizes, max_lin)
        else:
            np.testing.assert_equal(sizes[:tau_lin + 1], max_lin)
            block_size = tau_lin // 2
            i = block_size + 1
            for _ in range(2):
                j = i + block_size
                k = j + block_size
                np.testing.assert_allclose(sizes[i:j] / sizes[j:k], 2, atol=.5)
                i = j

    def check_pickling(self, acc):
        corr = acc.result()
        lags = acc.lag_times()
        sizes = acc.sample_sizes()
        acc_unpickled = pickle.loads(pickle.dumps(acc))
        np.testing.assert_array_equal(acc_unpickled.result(), corr)
        np.testing.assert_array_equal(acc_unpickled.lag_times(), lags)
        np.testing.assert_array_equal(acc_unpickled.sample_sizes(), sizes)

    def test_square_distance_componentwise(self):
        s = self.system
        v = np.array([1, 2, 3])
        p = s.part.add(pos=(0, 0, 0), v=v)

        obs = espressomd.observables.ParticlePositions(ids=(p.id,))
        acc = espressomd.accumulators.Correlator(
            obs1=obs, tau_lin=10, tau_max=2, delta_N=1,
            corr_operation="square_distance_componentwise")

        s.integrator.run(1000)
        s.auto_update_accumulators.add(acc)
        s.integrator.run(1000)

        # here don't call acc.finalize()
        corr = acc.result()
        self.check_pickling(acc)
        tau = self.calc_tau(s.time_step, acc.tau_lin, corr.shape[0])
        np.testing.assert_array_almost_equal(acc.lag_times(), tau)
        self.check_sizes(acc, 1000)
        for i in range(corr.shape[0]):
            np.testing.assert_array_almost_equal(corr[i], [v**2 * tau[i]**2])

    def test_tensor_product(self):
        s = self.system
        v = np.array([1, 2, 3])
        p = s.part.add(pos=(0, 0, 0), v=v)

        obs = espressomd.observables.ParticleVelocities(ids=(p.id,))
        acc = espressomd.accumulators.Correlator(
            obs1=obs, tau_lin=12, tau_max=2, delta_N=1,
            corr_operation="tensor_product")

        s.auto_update_accumulators.add(acc)
        s.integrator.run(1000)

        acc.finalize()
        corr = acc.result()
        self.check_pickling(acc)
        tau = self.calc_tau(s.time_step, acc.tau_lin, corr.shape[0])
        np.testing.assert_array_almost_equal(acc.lag_times(), tau)
        corr_ref = np.kron(v, v).reshape((3, 3))
        self.check_sizes(acc, 1000)
        for i in range(corr.shape[0]):
            np.testing.assert_array_almost_equal(corr[i], corr_ref)

    def test_componentwise_product(self):
        s = self.system
        v = np.array([1, 2, 3])
        p = s.part.add(pos=(0, 0, 0), v=v)

        obs = espressomd.observables.ParticleVelocities(ids=(p.id,))
        acc = espressomd.accumulators.Correlator(
            obs1=obs, tau_lin=10, tau_max=2, delta_N=1,
            corr_operation="componentwise_product")

        s.auto_update_accumulators.add(acc)
        s.integrator.run(1000)

        acc.finalize()
        corr = acc.result()
        self.check_pickling(acc)
        tau = self.calc_tau(s.time_step, acc.tau_lin, corr.shape[0])
        np.testing.assert_array_almost_equal(acc.lag_times(), tau)
        self.check_sizes(acc, 1000)
        for i in range(corr.shape[0]):
            np.testing.assert_array_almost_equal(corr[i], [v**2])

    def test_scalar_product(self):
        s = self.system
        v = np.array([1, 2, 3])
        p = s.part.add(pos=(0, 0, 0), v=v)

        obs = espressomd.observables.ParticleVelocities(ids=(p.id,))
        acc = espressomd.accumulators.Correlator(
            obs1=obs, tau_lin=10, tau_max=2, delta_N=1,
            corr_operation="scalar_product")

        s.auto_update_accumulators.add(acc)
        s.integrator.run(1000)

        acc.finalize()
        corr = acc.result()
        self.check_pickling(acc)
        tau = self.calc_tau(s.time_step, acc.tau_lin, corr.shape[0])
        np.testing.assert_array_almost_equal(acc.lag_times(), tau)
        self.check_sizes(acc, 1000)
        for i in range(corr.shape[0]):
            np.testing.assert_array_almost_equal(corr[i], [np.sum(v**2)])

    def test_fcs(self):
        s = self.system
        v = np.array([1, 2, 3])
        p = s.part.add(pos=(0, 0, 0), v=v)

        w = np.array([3, 2, 1])
        obs = espressomd.observables.ParticlePositions(ids=(p.id,))
        acc = espressomd.accumulators.Correlator(
            obs1=obs, tau_lin=10, tau_max=9.9 * self.system.time_step,
            delta_N=1, corr_operation="fcs_acf", args=w)

        s.auto_update_accumulators.add(acc)
        s.integrator.run(1000)

        acc.finalize()
        corr = acc.result()
        self.check_pickling(acc)
        tau = self.calc_tau(s.time_step, acc.tau_lin, corr.shape[0])
        np.testing.assert_array_almost_equal(acc.lag_times(), tau)
        self.check_sizes(acc, 1000, linear=True)
        for i in range(corr.shape[0]):
            np.testing.assert_array_almost_equal(
                corr[i],
                [np.exp(-np.linalg.norm(v / w * tau[i])**2)], decimal=10)

        # check setter and getter
        np.testing.assert_array_almost_equal(np.copy(acc.args), w**2)
        w_squared = np.array([4, 5, 6])**2
        acc.args = w_squared
        np.testing.assert_array_almost_equal(np.copy(acc.args), w_squared)

    def test_correlator_interface(self):
        # test setters and getters
        obs = espressomd.observables.ParticleVelocities(ids=(123,))
        acc = espressomd.accumulators.Correlator(
            obs1=obs, tau_lin=10, tau_max=12.0, delta_N=1,
            corr_operation="scalar_product")
        # check tau_lin
        self.assertEqual(acc.tau_lin, 10)
        # check tau_max
        self.assertEqual(acc.tau_max, 12.)
        # check delta_N
        self.assertEqual(acc.delta_N, 1)
        acc.delta_N = 2
        self.assertEqual(acc.delta_N, 2)
        # check corr_operation
        self.assertEqual(acc.corr_operation, "scalar_product")
        # check linear tau correlator and multiple tau correlator
        dt = self.system.time_step
        for tau_lin in (10, 20):
            for delta_N in (1, 2, 10):
                tau_max = dt * delta_N * tau_lin
                # linear, multiple and default (=multiple) tau correlator
                acc_lin = espressomd.accumulators.Correlator(
                    obs1=obs, tau_lin=tau_lin, tau_max=0.99 * tau_max,
                    delta_N=delta_N, corr_operation="scalar_product")
                acc_mul = espressomd.accumulators.Correlator(
                    obs1=obs, tau_lin=tau_lin, tau_max=1.0 * tau_max,
                    delta_N=delta_N, corr_operation="scalar_product")
                acc_def = espressomd.accumulators.Correlator(
                    obs1=obs, tau_lin=1, tau_max=tau_max,
                    delta_N=delta_N, corr_operation="scalar_product")
                lin_tau = acc_lin.lag_times()
                mul_tau = acc_mul.lag_times()
                def_tau = acc_mul.lag_times()
                # check tau
                time_lin = dt * delta_N * np.arange(len(lin_tau))
                time_mul = self.calc_tau(dt, tau_lin, len(mul_tau), delta_N)
                np.testing.assert_array_almost_equal(lin_tau, time_lin)
                np.testing.assert_array_almost_equal(mul_tau, time_mul)
                np.testing.assert_array_almost_equal(def_tau, time_mul)
                self.assertEqual(acc_def.tau_lin, tau_lin)
                # check pickling
                self.check_pickling(acc_lin)
                self.check_pickling(acc_lin)
                self.check_pickling(acc_def)


if __name__ == "__main__":
    ut.main()
