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
Testmodule for the observable recorder.

"""
import unittest as ut
import numpy as np
import espressomd
from espressomd.observables import ParticlePositions
from espressomd.accumulators import TimeSeries

N_PART = 100


class TimeSeriesTest(ut.TestCase):

    """
    Test class for the observable time series.

    """

    def test_time_series(self):
        """Check that accumulator results are the same as the respective numpy result.

        """

        system = espressomd.System(box_l=3 * [1.])
        system.part.add(pos=np.random.random((N_PART, 3)))

        obs = ParticlePositions(ids=system.part[:].id)
        time_series = TimeSeries(obs=obs)

        positions = []
        for _ in range(10):
            pos = np.random.random((N_PART, 3))
            positions.append(pos)

            system.part[:].pos = pos
            time_series.update()

        for result, expected in zip(time_series.time_series(), positions):
            np.testing.assert_array_equal(
                result, expected)

        time_series.clear()
        self.assertEqual(len(time_series.time_series()), 0)


if __name__ == "__main__":
    ut.main()
