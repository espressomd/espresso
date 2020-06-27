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

    def calc_tau(self, time_step, tau_lin, length):
        tau = []
        for i in range(tau_lin):
            tau.append(i)
        factor = 1
        while len(tau) < length:
            p = tau[-1] + factor * 1
            for i in range(0, tau_lin, 2):
                tau.append(p + factor * i)
            factor *= 2
        return time_step * np.array(tau[:length])

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
        acc_unpickeled = pickle.loads(pickle.dumps(acc))
        np.testing.assert_array_equal(corr, acc_unpickeled.result())

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
        acc_unpickeled = pickle.loads(pickle.dumps(acc))
        np.testing.assert_array_equal(corr, acc_unpickeled.result())

        tau = self.calc_tau(s.time_step, acc.tau_lin, corr.shape[0])
        np.testing.assert_array_almost_equal(corr[:, 0], tau)
        for i in range(corr.shape[0]):
            np.testing.assert_array_almost_equal(corr[i, 2:], np.sum(v**2))


if __name__ == "__main__":
    ut.main()
