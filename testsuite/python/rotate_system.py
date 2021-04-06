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
Testmodule for System.rotate_system()

"""
import unittest as ut
import numpy as np
import espressomd


class RotateSystemTest(ut.TestCase):
    system = espressomd.System(box_l=3 * [10.])

    def tearDown(self):
        self.system.part.clear()

    def test_no_mass(self):
        system = self.system
        p0 = system.part.add(pos=[4, 4, 4])
        p1 = system.part.add(pos=[6, 6, 6])

        pi = np.pi
        system.rotate_system(phi=0., theta=0., alpha=pi / 2.)

        np.testing.assert_allclose(np.copy(p0.pos), [6, 4, 4])
        np.testing.assert_allclose(np.copy(p1.pos), [4, 6, 6])

        system.rotate_system(phi=0., theta=0., alpha=-pi / 2.)

        np.testing.assert_allclose(np.copy(p0.pos), [4, 4, 4])
        np.testing.assert_allclose(np.copy(p1.pos), [6, 6, 6])

        system.rotate_system(phi=pi / 2., theta=0., alpha=pi / 2.)

        np.testing.assert_allclose(np.copy(p0.pos), [6, 4, 4])
        np.testing.assert_allclose(np.copy(p1.pos), [4, 6, 6])

        system.rotate_system(phi=pi / 2., theta=0., alpha=-pi / 2.)

        np.testing.assert_allclose(np.copy(p0.pos), [4, 4, 4])
        np.testing.assert_allclose(np.copy(p1.pos), [6, 6, 6])

        system.rotate_system(phi=pi / 2., theta=pi / 2., alpha=pi / 2.)

        np.testing.assert_allclose(np.copy(p0.pos), [4, 4, 6])
        np.testing.assert_allclose(np.copy(p1.pos), [6, 6, 4])

        system.rotate_system(phi=pi / 2., theta=pi / 2., alpha=-pi / 2.)

        np.testing.assert_allclose(np.copy(p0.pos), [4, 4, 4])
        np.testing.assert_allclose(np.copy(p1.pos), [6, 6, 6])

        # Check that virtual sites do not influence the center of mass
        # calculation
        if espressomd.has_features("VIRTUAL_SITES"):
            p2 = system.part.add(pos=p1.pos, virtual=True)
            system.rotate_system(phi=pi / 2., theta=pi / 2., alpha=-pi / 2.)
            np.testing.assert_allclose(np.copy(p0.pos), [6, 4, 4])
            np.testing.assert_allclose(np.copy(p1.pos), [4, 6, 6])
            np.testing.assert_allclose(np.copy(p2.pos), [4, 6, 6])


if __name__ == "__main__":
    ut.main()
