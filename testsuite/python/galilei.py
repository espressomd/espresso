from __future__ import print_function

import unittest as ut
import numpy as np

import espressomd
from espressomd.galilei import GalileiTransform

BOX_L = np.array([10, 20, 30])
N_PART = 500


class Galilei(ut.TestCase):
    system = espressomd.System(box_l=BOX_L)

    def setUp(self):
        self.system.part.add(pos=BOX_L * np.random.random((N_PART, 3)),
                             v=-5. + 10. * np.random.random((N_PART, 3)),
                             f=np.random.random((N_PART, 3)))
        if espressomd.has_features("MASS"):
            self.system.part[:].mass = 42. * np.random.random((N_PART,))

    def tearDown(self):
        self.system.part.clear()

    def test_kill_particle_motion(self):
        g = GalileiTransform()
        g.kill_particle_motion()

        np.testing.assert_array_equal(
            np.copy(self.system.part[:].v), np.zeros((N_PART, 3)))

    def test_kill_particle_forces(self):
        g = GalileiTransform()
        g.kill_particle_forces()

        np.testing.assert_array_equal(
            np.copy(self.system.part[:].f), np.zeros((N_PART, 3)))

    def test_cms(self):
        parts = self.system.part[:]
        g = GalileiTransform()

        total_mass = np.sum(parts.mass)
        com = np.sum(
            np.multiply(parts.mass.reshape((N_PART, 1)), parts.pos), axis=0) / total_mass

        np.testing.assert_allclose(np.copy(g.system_CMS()), com)

    def test_cms_velocity(self):
        parts = self.system.part[:]
        g = GalileiTransform()
        total_mass = np.sum(parts.mass)
        com_v = np.sum(
            np.multiply(parts.mass.reshape((N_PART, 1)), parts.v), axis=0) / total_mass

        np.testing.assert_allclose(np.copy(g.system_CMS_velocity()), com_v)

    def test_galilei_transform(self):
        g = GalileiTransform()
        g.galilei_transform()

        np.testing.assert_allclose(
            np.copy(g.system_CMS_velocity()), np.zeros((3,)), atol=1e-15)

if __name__ == "__main__":
    ut.main()
