# Copyright (C) 2019 The ESPResSo project
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

import espressomd
import unittest_decorators as utx

N_PART = 10
VELOCITY = np.array([1.0, 2.0, 3.0])
MASS = 2.1


@utx.skipIfMissingFeatures("MASS")
class LinearMomentumTest(ut.TestCase):
    system = None

    @classmethod
    def setUpClass(cls):
        cls.system = espressomd.System(box_l=[10.0] * 3)

    def test(self):
        self.system.part.add(pos=np.random.random((N_PART, 3)), v=np.ones(
            (N_PART, 3)) * VELOCITY, mass=np.ones((N_PART)) * MASS)
        linear_momentum = self.system.analysis.linear_momentum()
        np.testing.assert_allclose(linear_momentum, N_PART * MASS * VELOCITY)


if __name__ == "__main__":
    ut.main()
