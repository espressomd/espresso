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
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])

    def setUp(self):
        self.system.box_l = [10, 10, 10]
        self.system.part.add(id=0, pos=(0.1, 0.1, 0.1), dip=(1.3, 2.1, -6))
        self.system.part.add(id=1, pos=(0, 0, 0), dip=(7.3, 6.1, -4))

    def tearDown(self):
        self.system.part.clear()
        self.system.actors.clear()

    if espressomd.has_features(["DP3M"]):
        test_DP3M = generate_test_for_class(
            system, magnetostatics.DipolarP3M,
            dict(prefactor=1.0, epsilon=0.0, inter=1000,
                 mesh_off=[0.5, 0.5, 0.5], r_cut=2.4, mesh=[8, 8, 8],
                 cao=1, alpha=12, accuracy=0.01, tune=False))

    if espressomd.has_features(["DIPOLAR_DIRECT_SUM"]):
        test_DdsCpu = generate_test_for_class(
            system, magnetostatics.DipolarDirectSumCpu, dict(prefactor=3.4))
        if espressomd.has_features("EXPERIMENTAL_FEATURES"):
            test_DdsRCpu = generate_test_for_class(
                system, magnetostatics.DipolarDirectSumWithReplicaCpu,
                dict(prefactor=3.4, n_replica=2))

    @utx.skipIfMissingFeatures(["DP3M"])
    def test_dp3m(self):
        self.system.time_step = 0.01
        prefactor = 1.1
        self.system.part[0].pos = [1.0, 2.0, 2.0]
        self.system.part[1].pos = [3.0, 2.0, 2.0]
        dp3m_energy = 1.6733349639532644
        dp3m_force = np.array([-3.54175042, -4.6761059, 9.96632774])
        dp3m = espressomd.magnetostatics.DipolarP3M(
            prefactor=prefactor,
            accuracy=9.995178689932661e-07,
            mesh=[49, 49, 49],
            mesh_off=[0.5, 0.5, 0.5],
            cao=7,
            r_cut=4.739799499511719,
            alpha=0.9056147262573242,
            alpha_L=9.056147262573242,
            r_cut_iL=0.4739799499511719,
            cao_cut=3 * [0.7142857142857142],
            a=3 * [0.2040816326530612],
            additional_mesh=[0., 0., 0.],
            tune=False)
        self.system.actors.add(dp3m)
        self.assertAlmostEqual(self.system.analysis.energy()['dipolar'],
                               dp3m_energy, places=5)
        # need to update forces
        self.system.integrator.run(0)
        np.testing.assert_allclose(np.copy(self.system.part[0].f),
                                   dp3m_force, atol=1E-5)
        np.testing.assert_allclose(np.copy(self.system.part[1].f),
                                   -dp3m_force, atol=1E-5)

    @utx.skipIfMissingFeatures(["DP3M"])
    def test_dp3m_non_metallic(self):
        self.system.time_step = 0.01
        prefactor = 1.1
        self.system.part[0].pos = [1.0, 2.0, 2.0]
        self.system.part[1].pos = [3.0, 2.0, 2.0]
        for epsilon_power in range(-4, 5):
            epsilon = 10**epsilon_power
            dp3m_energy = 1.66706 / (1 + 2 * epsilon) + 1.673333
            dp3m = espressomd.magnetostatics.DipolarP3M(
                prefactor=prefactor,
                accuracy=9.995178689932661e-07,
                mesh=[49, 49, 49],
                mesh_off=[0.5, 0.5, 0.5],
                cao=7,
                epsilon=epsilon,
                r_cut=4.739799499511719,
                alpha=0.9056147262573242,
                alpha_L=9.056147262573242,
                r_cut_iL=0.4739799499511719,
                cao_cut=3 * [0.7142857142857142],
                a=3 * [0.2040816326530612],
                additional_mesh=[0., 0., 0.],
                tune=False)
            self.system.actors.add(dp3m)
            self.assertAlmostEqual(self.system.analysis.energy()['dipolar'],
                                   dp3m_energy, places=3)
            self.system.actors.remove(dp3m)


if __name__ == "__main__":
    ut.main()
