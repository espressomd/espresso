#
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
import espressomd.galilei


class Galilei(ut.TestCase):
    system = espressomd.System(box_l=[10, 20, 30])

    def setUp(self):
        N_PART = 500
        self.partcls = self.system.part.add(
            pos=self.system.box_l * np.random.random((N_PART, 3)),
            v=-5. + 10. * np.random.random((N_PART, 3)),
            f=np.random.random((N_PART, 3)))
        if espressomd.has_features("MASS"):
            self.partcls.mass = 42. * np.random.random((N_PART,))
        if espressomd.has_features("ROTATION"):
            self.partcls.omega_lab = -2. - np.random.random((N_PART, 3))
            self.partcls.torque_lab = - \
                2. - np.random.random((N_PART, 3))

    def tearDown(self):
        self.system.part.clear()

    def test_kill_particle_motion(self):
        g = espressomd.galilei.GalileiTransform()
        g.kill_particle_motion()

        np.testing.assert_array_equal(np.copy(self.partcls.v), 0)

        if espressomd.has_features("ROTATION"):
            np.testing.assert_array_less(
                np.copy(self.partcls.omega_lab), 0)

            g.kill_particle_motion(rotation=True)

            np.testing.assert_array_equal(
                np.copy(self.partcls.omega_lab), 0)

    def test_kill_particle_forces(self):
        g = espressomd.galilei.GalileiTransform()
        g.kill_particle_forces()

        np.testing.assert_array_equal(np.copy(self.partcls.f), 0)

        if espressomd.has_features("ROTATION"):
            np.testing.assert_array_less(
                np.copy(self.partcls.torque_lab), 0)

            g.kill_particle_forces(torque=True)

            np.testing.assert_array_equal(
                np.copy(self.partcls.torque_lab), 0)

    def test_cms(self):
        parts = self.partcls
        g = espressomd.galilei.GalileiTransform()

        total_mass = np.sum(parts.mass)
        com = np.sum(
            np.multiply(parts.mass.reshape((-1, 1)), parts.pos), axis=0) / total_mass

        np.testing.assert_allclose(np.copy(g.system_CMS()), com)

    def test_cms_velocity(self):
        parts = self.partcls
        g = espressomd.galilei.GalileiTransform()
        total_mass = np.sum(parts.mass)
        com_v = np.sum(
            np.multiply(parts.mass.reshape((-1, 1)), parts.v), axis=0) / total_mass

        np.testing.assert_allclose(np.copy(g.system_CMS_velocity()), com_v)

    def test_galilei_transform(self):
        g = espressomd.galilei.GalileiTransform()
        g.galilei_transform()

        np.testing.assert_allclose(
            np.copy(g.system_CMS_velocity()), np.zeros((3,)), atol=1e-14)


if __name__ == "__main__":
    ut.main()
