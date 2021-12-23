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

import unittest as ut

import numpy as np
import pickle

import espressomd
import espressomd.observables
import espressomd.accumulators

N_PART = 4


class TimeSeriesTest(ut.TestCase):

    """
    Test class for the TimeSeries accumulator.

    """
    system = espressomd.System(box_l=[10.0] * 3)
    system.cell_system.skin = 0.4
    system.time_step = 0.01

    def setUp(self):
        np.random.seed(seed=42)

    def tearDown(self):
        self.system.part.clear()
        self.system.auto_update_accumulators.clear()

    def test_accumulator(self):
        """Check that accumulator results are the same as the respective numpy result.

        """
        system = self.system
        system.part.add(pos=np.zeros((N_PART, 3)))
        obs = espressomd.observables.ParticlePositions(ids=range(N_PART))
        acc = espressomd.accumulators.TimeSeries(obs=obs)
        system.auto_update_accumulators.add(acc)
        positions = np.copy(system.box_l) * np.random.random((10, N_PART, 3))

        for pos in positions:
            system.part.all().pos = pos
            acc.update()

        self.assertEqual(acc.get_params()['obs'], obs)

        np.testing.assert_array_equal(acc.time_series(), positions)

        # Check pickling
        acc_unpickled = pickle.loads(pickle.dumps(acc))
        np.testing.assert_array_equal(acc_unpickled.time_series(), positions)

        acc.clear()
        self.assertEqual(len(acc.time_series()), 0)


if __name__ == "__main__":
    ut.main()
