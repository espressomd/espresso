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

    def test_correlator_compression(self):
        p = self.system.part.add(pos=(0, 0, 0))
        obs = espressomd.observables.ParticleVelocities(ids=(0,))
        v1 = 3.
        v2 = 5.
        compressed_ref = {
            "discard1": v2**2, "linear": ((v1 + v2) / 2.)**2,
            "discard2": v1**2, "uncompressed": ((v1**2 + v2**2) / 2., v1 * v2)}
        for tau_lin in [8, 10, 12]:
            # set up accumulators
            accumulators = {}
            for compression in ("linear", "discard1", "discard2"):
                acc = espressomd.accumulators.Correlator(
                    obs1=obs, tau_lin=tau_lin, tau_max=2, delta_N=1,
                    compress1=compression, compress2=compression,
                    corr_operation="scalar_product")
                accumulators[compression] = acc
                self.system.auto_update_accumulators.add(acc)
            # record oscillating data with frequency of 1/time_step
            for i in range(1000):
                p.v = (v1 if (i % 2 == 0) else v2, 0, 0)
                self.system.integrator.run(1)
            # check compression algorithms: the first tau_lin values
            # are always uncompressed, the rest is compressed; the
            # oscillation frequency is such that 'discard*' methods
            # yield the corresponding constant 'v*' while 'linear'
            # yields the average of 'v1' and 'v2'
            uncompressed = np.repeat(
                [compressed_ref["uncompressed"]], tau_lin, axis=0).flatten()
            for compression in ("linear", "discard1", "discard2"):
                accumulators[compression].finalize()
                corr = accumulators[compression].result().flatten()
                corr_ref = np.repeat(compressed_ref[compression], corr.shape)
                corr_ref[:tau_lin + 1] = uncompressed[:tau_lin + 1]
                np.testing.assert_array_equal(corr, corr_ref)
            self.system.auto_update_accumulators.clear()

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

    def test_correlator_exceptions(self):
        self.system.part.add(pos=2 * [(0, 0, 0)])
        obs = espressomd.observables.ParticleVelocities(ids=(0,))

        def create_accumulator(obs1=obs, **kwargs):
            valid_kwargs = {'obs1': obs1, 'tau_lin': 10, 'tau_max': 10.,
                            'delta_N': 1, 'corr_operation': "scalar_product"}
            valid_kwargs.update(kwargs)
            return espressomd.accumulators.Correlator(**valid_kwargs)

        # check finalize method
        acc = create_accumulator()
        acc.finalize()
        with self.assertRaisesRegex(RuntimeError, r"Correlator::finalize\(\) can only be called once"):
            acc.finalize()
        with self.assertRaisesRegex(RuntimeError, r"No data can be added after finalize\(\) was called."):
            acc.update()

        # check general arguments and input data
        with self.assertRaisesRegex(RuntimeError, "tau_lin must be >= 2"):
            create_accumulator(tau_lin=0)
        with self.assertRaisesRegex(RuntimeError, "tau_lin must be divisible by 2"):
            create_accumulator(tau_lin=3)
        with self.assertRaisesRegex(RuntimeError, "tau_max must be >= delta_t"):
            create_accumulator(delta_N=2 * int(10. / self.system.time_step))
        with self.assertRaisesRegex(ValueError, "correlation operation 'unknown' not implemented"):
            create_accumulator(corr_operation="unknown")
        with self.assertRaisesRegex(ValueError, "unknown compression method 'unknown1' for first observable"):
            create_accumulator(compress1="unknown1")
        with self.assertRaisesRegex(ValueError, "unknown compression method 'unknown2' for second observable"):
            create_accumulator(compress2="unknown2")
        with self.assertRaisesRegex(RuntimeError, "dimension of first observable has to be >= 1"):
            create_accumulator(
                obs1=espressomd.observables.ParticleVelocities(ids=()))
        with self.assertRaisesRegex(RuntimeError, "dimension of second observable has to be >= 1"):
            create_accumulator(
                obs2=espressomd.observables.ParticleVelocities(ids=()))

        # check FCS-specific arguments and input data
        with self.assertRaisesRegex(RuntimeError, "missing parameter for fcs_acf: w_x w_y w_z"):
            create_accumulator(corr_operation="fcs_acf", args=[1, 1, 0])
        with self.assertRaisesRegex(RuntimeError, "dimA must be divisible by 3 for fcs_acf"):
            create_accumulator(corr_operation="fcs_acf", args=[1, 1, 1],
                               obs1=espressomd.observables.Energy())
        with self.assertRaisesRegex(RuntimeError, "the last dimension of dimA must be 3 for fcs_acf"):
            obs_dens = espressomd.observables.DensityProfile(
                ids=(0,), n_x_bins=3, n_y_bins=3, n_z_bins=1, min_x=0.,
                min_y=0., min_z=0., max_x=1., max_y=1., max_z=1.)
            create_accumulator(corr_operation="fcs_acf", args=[1, 1, 1],
                               obs1=obs_dens)

        # check correlation errors
        obs2 = espressomd.observables.ParticleVelocities(ids=(0, 1))
        with self.assertRaisesRegex(RuntimeError, "Error in scalar product: The vector sizes do not match"):
            acc = create_accumulator(obs2=obs2)
            acc.update()
        with self.assertRaisesRegex(RuntimeError, "Error in componentwise product: The vector sizes do not match"):
            acc = create_accumulator(
                obs2=obs2, corr_operation="componentwise_product")
            acc.update()
        with self.assertRaisesRegex(RuntimeError, "Error in square distance componentwise: The vector sizes do not match"):
            acc = create_accumulator(
                obs2=obs2, corr_operation="square_distance_componentwise")
            acc.update()
        with self.assertRaisesRegex(RuntimeError, "Error in fcs_acf: The vector sizes do not match"):
            acc = create_accumulator(
                obs2=obs2, corr_operation="fcs_acf", args=[1, 1, 1])
            acc.update()


if __name__ == "__main__":
    ut.main()
