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
    # Handle for espresso system
    system = espressomd.System(box_l=[10, 10, 10])
    system.cell_system.skin = 0.4
    system.time_step = 0.01

    def tearDown(self):
        self.system.part.clear()

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

    def test_square_distance_componentwise(self):
        s = self.system
        s.part.add(id=0, pos=(0, 0, 0), v=(1, 2, 3))

        obs = espressomd.observables.ParticlePositions(ids=(0,))
        acc = espressomd.accumulators.Correlator(
            obs1=obs, tau_lin=10, tau_max=10.0, delta_N=1,
            corr_operation="square_distance_componentwise")

        s.integrator.run(1000)
        s.auto_update_accumulators.add(acc)
        s.integrator.run(20000)

        corr = acc.result()

        # Check pickling
        acc_unpickled = pickle.loads(pickle.dumps(acc))
        np.testing.assert_array_equal(corr, acc_unpickled.result())

        tau = self.calc_tau(s.time_step, acc.tau_lin, corr.shape[0])
        np.testing.assert_array_almost_equal(corr[:, 0], tau)
        for i in range(corr.shape[0]):
            t = corr[i, 0]
            self.assertAlmostEqual(corr[i, 2], t * t, places=3)
            self.assertAlmostEqual(corr[i, 3], 4 * t * t, places=3)
            self.assertAlmostEqual(corr[i, 4], 9 * t * t, places=3)

    def test_tensor_product(self):
        s = self.system
        v = np.array([1, 2, 3])
        s.part.add(id=0, pos=(0, 0, 0), v=v)

        obs = espressomd.observables.ParticleVelocities(ids=(0,))
        acc = espressomd.accumulators.Correlator(
            obs1=obs, tau_lin=10, tau_max=10.0, delta_N=1,
            corr_operation="tensor_product")

        s.integrator.run(1000)
        s.auto_update_accumulators.add(acc)
        s.integrator.run(20000)

        corr = acc.result()

        # Check pickling
        acc_unpickled = pickle.loads(pickle.dumps(acc))
        np.testing.assert_array_equal(corr, acc_unpickled.result())

        tau = self.calc_tau(s.time_step, acc.tau_lin, corr.shape[0])
        np.testing.assert_array_almost_equal(corr[:, 0], tau)
        for i in range(corr.shape[0]):
            np.testing.assert_array_almost_equal(corr[i, 2:], np.kron(v, v))

    def test_componentwise_product(self):
        s = self.system
        v = np.array([1, 2, 3])
        s.part.add(id=0, pos=(0, 0, 0), v=v)

        obs = espressomd.observables.ParticleVelocities(ids=(0,))
        acc = espressomd.accumulators.Correlator(
            obs1=obs, tau_lin=10, tau_max=10.0, delta_N=1,
            corr_operation="componentwise_product")

        s.integrator.run(1000)
        s.auto_update_accumulators.add(acc)
        s.integrator.run(20000)

        corr = acc.result()

        # Check pickling
        acc_unpickled = pickle.loads(pickle.dumps(acc))
        np.testing.assert_array_equal(corr, acc_unpickled.result())

        tau = self.calc_tau(s.time_step, acc.tau_lin, corr.shape[0])
        np.testing.assert_array_almost_equal(corr[:, 0], tau)
        for i in range(corr.shape[0]):
            np.testing.assert_array_almost_equal(corr[i, 2:], v**2)

    def test_scalar_product(self):
        s = self.system
        v = np.array([1, 2, 3])
        s.part.add(id=0, pos=(0, 0, 0), v=v)

        obs = espressomd.observables.ParticleVelocities(ids=(0,))
        acc = espressomd.accumulators.Correlator(
            obs1=obs, tau_lin=10, tau_max=10.0, delta_N=1,
            corr_operation="scalar_product")

        s.integrator.run(1000)
        s.auto_update_accumulators.add(acc)
        s.integrator.run(20000)

        corr = acc.result()

        # Check pickling
        acc_unpickled = pickle.loads(pickle.dumps(acc))
        np.testing.assert_array_equal(corr, acc_unpickled.result())

        tau = self.calc_tau(s.time_step, acc.tau_lin, corr.shape[0])
        np.testing.assert_array_almost_equal(corr[:, 0], tau)
        for i in range(corr.shape[0]):
            np.testing.assert_array_almost_equal(corr[i, 2:], np.sum(v**2))

    def test_fcs(self):
        s = self.system
        v = np.array([1, 2, 3])
        s.part.add(id=0, pos=(0, 0, 0), v=v)

        w = np.array([3, 2, 1])
        obs = espressomd.observables.ParticlePositions(ids=(0,))
        acc = espressomd.accumulators.Correlator(
            obs1=obs, tau_lin=10, tau_max=2.0, delta_N=1,
            corr_operation="fcs_acf", args=w)

        s.integrator.run(100)
        s.auto_update_accumulators.add(acc)
        s.integrator.run(2000)

        corr = acc.result()

        # Check pickling
        acc_unpickled = pickle.loads(pickle.dumps(acc))
        np.testing.assert_array_equal(corr, acc_unpickled.result())

        tau = self.calc_tau(s.time_step, acc.tau_lin, corr.shape[0])
        np.testing.assert_array_almost_equal(corr[:, 0], tau)
        for i in range(corr.shape[0]):
            np.testing.assert_array_almost_equal(
                corr[i, 2:],
                np.exp(-np.linalg.norm(v / w * tau[i])**2), decimal=10)

        # check setter and getter
        np.testing.assert_array_almost_equal(np.copy(acc.args), w**2)
        w_squared = np.array([4, 5, 6])**2
        acc.args = w_squared
        np.testing.assert_array_almost_equal(np.copy(acc.args), w_squared)

    def test_correlator_interface(self):
        # test setters and getters
        obs = espressomd.observables.ParticleVelocities(ids=(0,))
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
                corr_lin = acc_lin.result()
                corr_mul = acc_mul.result()
                corr_def = acc_mul.result()
                # check tau
                time_lin = dt * delta_N * np.arange(len(corr_lin))
                time_mul = self.calc_tau(dt, tau_lin, len(corr_mul), delta_N)
                np.testing.assert_array_almost_equal(corr_lin[:, 0], time_lin)
                np.testing.assert_array_almost_equal(corr_mul[:, 0], time_mul)
                np.testing.assert_array_almost_equal(corr_def[:, 0], time_mul)
                self.assertEqual(acc_def.tau_lin, tau_lin)
                # check pickling
                corr_lin_upkl = pickle.loads(pickle.dumps(acc_lin))
                corr_mul_upkl = pickle.loads(pickle.dumps(acc_mul))
                np.testing.assert_array_equal(corr_lin, corr_lin_upkl.result())
                np.testing.assert_array_equal(corr_mul, corr_mul_upkl.result())


if __name__ == "__main__":
    ut.main()
