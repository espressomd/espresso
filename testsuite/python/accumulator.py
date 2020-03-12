#
# Copyright (C) 2017-2019 The ESPResSo project
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

"""
Testmodule for the observable accumulator.

"""
import sys
import unittest as ut
import numpy as np
import espressomd
import espressomd.observables
import espressomd.accumulators


N_PART = 4


class AccumulatorTest(ut.TestCase):

    """
    Test class for the observable accumulator.

    """

    def setUp(self):
        np.random.seed(seed=162)
        self.system = espressomd.System(box_l=[10.0] * 3)
        self.system.cell_system.skin = 0.4
        self.system.time_step = 0.01
        self.system.part.add(pos=np.zeros((N_PART, 3)))
        self.system.integrator.run(steps=0)
        self.pos_obs = espressomd.observables.ParticlePositions(
            ids=range(N_PART))
        self.pos_obs_acc = espressomd.accumulators.MeanVarianceCalculator(
            obs=self.pos_obs)
        self.system.auto_update_accumulators.add(self.pos_obs_acc)
        self.positions = np.copy(self.system.box_l) * \
            np.random.random((10, N_PART, 3))

    def test_accumulator(self):
        """Check that accumulator results are the same as the respective numpy result.

        """
        for pos in self.positions:
            self.system.part[:].pos = pos
            self.system.integrator.run(1)
        self.assertEqual(self.pos_obs, self.pos_obs_acc.get_params()['obs'])
        acc_mean = self.pos_obs_acc.get_mean()
        np.testing.assert_allclose(
            acc_mean,
            np.mean(self.positions, axis=0), atol=1e-4)
        np.testing.assert_allclose(
            self.pos_obs_acc.get_variance(),
            np.var(self.positions, axis=0, ddof=1), atol=1e-4)


if __name__ == "__main__":
    suite = ut.TestSuite()
    suite.addTests(ut.TestLoader().loadTestsFromTestCase(AccumulatorTest))
    result = ut.TextTestRunner(verbosity=4).run(suite)
    sys.exit(not result.wasSuccessful())
