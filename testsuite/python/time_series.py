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
Testmodule for the time series accumulator.

"""
import unittest as ut
import numpy as np
import pickle

import espressomd
import espressomd.observables
import espressomd.accumulators

N_PART = 100


class TimeSeriesTest(ut.TestCase):

    def test_time_series(self):
        """Check that accumulator results are the same as the respective numpy result.

        """

        system = espressomd.System(box_l=3 * [1.])
        system.part.add(pos=np.random.random((N_PART, 3)))

        obs = espressomd.observables.ParticlePositions(ids=system.part[:].id)
        acc = espressomd.accumulators.TimeSeries(obs=obs)

        positions = []
        for _ in range(10):
            pos = np.random.random((N_PART, 3))
            positions.append(pos)

            system.part[:].pos = pos
            acc.update()

        time_series = acc.time_series()

        # Check pickling
        acc_unpickled = pickle.loads(pickle.dumps(acc))
        np.testing.assert_array_equal(time_series, acc_unpickled.time_series())

        for result, expected in zip(time_series, positions):
            np.testing.assert_array_equal(
                np.array(result).reshape((N_PART, 3)), expected)

        acc.clear()
        self.assertEqual(len(acc.time_series()), 0)


if __name__ == "__main__":
    ut.main()
