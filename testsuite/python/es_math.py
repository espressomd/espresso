# Copyright (C) 2010-2019 The ESPResSo project
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
import numpy as np
import unittest as ut
import espressomd.math


class TestMath(ut.TestCase):

    def check_orthonormality(self, vec1, vec2):
        self.assertAlmostEqual(np.linalg.norm(vec1), 1.)
        self.assertAlmostEqual(np.linalg.norm(vec2), 1.)
        self.assertAlmostEqual(np.dot(vec1, vec2), 0)

    def test_cylindrical_transformation_parameters(self):
        """ Test for the varous constructors of CylindricalTransformationParameters """

        ctp_default = espressomd.math.CylindricalTransformationParameters()
        self.check_orthonormality(ctp_default.axis, ctp_default.orientation)

        axis = np.array([-17, 0.1, np.pi])
        axis /= np.linalg.norm(axis)
        ctp_auto_orientation = espressomd.math.CylindricalTransformationParameters(
            center=3 * [42], axis=axis)
        self.check_orthonormality(
            ctp_auto_orientation.axis,
            ctp_auto_orientation.orientation)

        ctp_full = espressomd.math.CylindricalTransformationParameters(
            center=3 * [42], axis=[0, 1, 0], orientation=[1, 0, 0])
        self.check_orthonormality(ctp_full.axis, ctp_full.orientation)

        with self.assertRaises(Exception):
            ctp_only_center = espressomd.math.CylindricalTransformationParameters(
                center=3 * [42])
            ctp_only_center.axis = 3 * [3]


if __name__ == "__main__":
    ut.main()
