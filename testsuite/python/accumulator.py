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
    system = espressomd.System(box_l=[10.0] * 3)
    system.cell_system.skin = 0.4
    system.time_step = 0.01

    def setUp(self):
        np.random.seed(seed=42)

    def test_accumulator(self):
        """Check that accumulator results are the same as the respective numpy result.

        """
        system = self.system
        system.part.add(pos=np.zeros((N_PART, 3)))
        system.integrator.run(steps=0)
        pos_obs = espressomd.observables.ParticlePositions(
            ids=range(N_PART))
        pos_obs_acc = espressomd.accumulators.MeanVarianceCalculator(
            obs=pos_obs)
        system.auto_update_accumulators.add(pos_obs_acc)
        positions = np.copy(system.box_l) * \
            np.random.random((10, N_PART, 3))
        for pos in positions:
            system.part[:].pos = pos
            system.integrator.run(1)
        self.assertEqual(pos_obs, pos_obs_acc.get_params()['obs'])
        np.testing.assert_allclose(
            pos_obs_acc.mean(),
            np.mean(positions, axis=0), atol=1e-4)
        variance = np.var(positions, axis=0, ddof=1)
        np.testing.assert_allclose(
            pos_obs_acc.variance(),
            variance, atol=1e-4)
        np.testing.assert_allclose(
            pos_obs_acc.std_error(),
            np.sqrt(variance / len(positions)), atol=1e-4)


if __name__ == "__main__":
    ut.main()
