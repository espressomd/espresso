# Copyright (C) 2010-2018 The ESPResSo project
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
import sys
import numpy as np
import unittest as ut

from espressomd import utils


class TestUnravelIndex(ut.TestCase):

    def test_against_numpy(self):
        xdim = (0, 10)
        ydim = (0, 200)
        zdim = (0, 150)
        edim = (2, 244)
        index = 4127
        numpy_result = np.unravel_index(
            index,
            (xdim[1] - xdim[0],
             ydim[1] - ydim[0],
                zdim[1] - zdim[0],
                edim[1] - edim[0]))
        utils_result = np.array(utils.get_unravelled_index(
            [xdim[1] - xdim[0], ydim[1] - ydim[0], zdim[1] - zdim[0], edim[1] - edim[0]], 4, index))
        np.testing.assert_array_equal(numpy_result, utils_result)


if __name__ == "__main__":
    suite = ut.TestSuite()
    suite.addTests(ut.TestLoader().loadTestsFromTestCase(
        TestUnravelIndex))
    result = ut.TextTestRunner(verbosity=4).run(suite)
    sys.exit(not result.wasSuccessful())
