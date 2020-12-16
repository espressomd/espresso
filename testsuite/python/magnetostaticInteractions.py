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
import unittest_decorators as utx
import numpy as np

import espressomd
from espressomd import magnetostatics
from tests_common import generate_test_for_class


class MagnetostaticsInteractionsTests(ut.TestCase):
    # Handle to espresso system
    system = espressomd.System(box_l=[10., 10., 10.])

    def setUp(self):
        self.system.box_l = [10., 10., 10.]
        self.system.part.add(id=0, pos=(0.1, 0.1, 0.1), dip=(1.3, 2.1, -6))
        self.system.part.add(id=1, pos=(0, 0, 0), dip=(7.3, 6.1, -4))

    def tearDown(self):
        self.system.part.clear()
        self.system.actors.clear()

    if espressomd.has_features(["DP3M"]):
        test_DP3M = generate_test_for_class(
            system, magnetostatics.DipolarP3M,
            dict(prefactor=1., epsilon=0., mesh_off=[0.5, 0.5, 0.5], r_cut=2.4,
                 cao=1, mesh=[8, 8, 8], alpha=12, accuracy=0.01, tune=False))

    test_DdsCpu = generate_test_for_class(
        system, magnetostatics.DipolarDirectSumCpu, dict(prefactor=3.4))

    test_DdsRCpu = generate_test_for_class(
        system, magnetostatics.DipolarDirectSumWithReplicaCpu,
        dict(prefactor=3.4, n_replica=2))

    def ref_values(self, epsilon=np.inf):
        x = 1. / (1 + 2 * epsilon)
        dp3m_energy = 1.66706 * x + 1.673333
        dp3m_torque1 = np.array([-0.5706503 * x + 2.561371,
                                 -0.1812375 * x + 10.394144,
                                 -0.2976916 * x + 9.965342])
        dp3m_torque2 = np.array([+0.3362938 * x + 1.854679,
                                 -0.2269749 * x - 3.638175,
                                 +0.5315054 * x + 8.487292])
        dp3m_force = np.array([-3.54175042, -4.6761059, 9.96632774])
        alpha, r_cut, mesh, cao = (9.056147262573242, 4.739799499511719, 49, 7)
        dp3m_params = {'prefactor': 1.1, 'accuracy': 9.995178689932661e-07,
                       'mesh': mesh, 'mesh_off': [0.5, 0.5, 0.5],
                       'cao': cao, 'additional_mesh': [0.0, 0.0, 0.0],
                       'alpha': alpha / 10, 'alpha_L': alpha, 'r_cut': r_cut,
                       'r_cut_iL': r_cut / self.system.box_l[0],
                       'cao_cut': 3 * [self.system.box_l[0] / mesh / 2 * cao],
                       'a': 3 * [self.system.box_l[0] / mesh]}
        return dp3m_params, dp3m_energy, dp3m_force, dp3m_torque1, dp3m_torque2

    @utx.skipIfMissingFeatures(["DP3M"])
    def test_dp3m(self):
        self.system.time_step = 0.01
        self.system.part[0].pos = [1.0, 2.0, 2.0]
        self.system.part[1].pos = [3.0, 2.0, 2.0]
        dp3m_params, dp3m_energy, dp3m_force, dp3m_torque1, dp3m_torque2 = self.ref_values()
        dp3m = espressomd.magnetostatics.DipolarP3M(tune=False, **dp3m_params)
        self.system.actors.add(dp3m)
        self.assertAlmostEqual(self.system.analysis.energy()['dipolar'],
                               dp3m_energy, places=5)
        # update forces and torques
        self.system.integrator.run(0)
        np.testing.assert_allclose(np.copy(self.system.part[0].f),
                                   dp3m_force, atol=1E-5)
        np.testing.assert_allclose(np.copy(self.system.part[1].f),
                                   -dp3m_force, atol=1E-5)
        np.testing.assert_allclose(np.copy(self.system.part[0].torque_lab),
                                   dp3m_torque1, atol=1E-5)
        np.testing.assert_allclose(np.copy(self.system.part[1].torque_lab),
                                   dp3m_torque2, atol=1E-5)

    @utx.skipIfMissingFeatures(["DP3M"])
    def test_dp3m_non_metallic(self):
        self.system.time_step = 0.01
        self.system.part[0].pos = [1.0, 2.0, 2.0]
        self.system.part[1].pos = [3.0, 2.0, 2.0]
        for epsilon_power in range(-4, 5):
            epsilon = 10**epsilon_power
            dp3m_params, dp3m_energy, dp3m_force, dp3m_torque1, dp3m_torque2 = self.ref_values(
                epsilon)
            dp3m = espressomd.magnetostatics.DipolarP3M(
                tune=False, epsilon=epsilon, **dp3m_params)
            self.system.actors.add(dp3m)
            self.assertAlmostEqual(self.system.analysis.energy()['dipolar'],
                                   dp3m_energy, places=5)
            # update forces and torques
            self.system.integrator.run(0)
            np.testing.assert_allclose(np.copy(self.system.part[0].f),
                                       dp3m_force, atol=1E-5)
            np.testing.assert_allclose(np.copy(self.system.part[1].f),
                                       -dp3m_force, atol=1E-5)
            np.testing.assert_allclose(np.copy(self.system.part[0].torque_lab),
                                       dp3m_torque1, atol=1E-5)
            np.testing.assert_allclose(np.copy(self.system.part[1].torque_lab),
                                       dp3m_torque2, atol=1E-5)
            self.system.actors.remove(dp3m)


if __name__ == "__main__":
    ut.main()
