#
# Copyright (C) 2013-2019 The ESPResSo project
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
import unittest_decorators as utx
import espressomd
import numpy as np
import itertools


@utx.skipIfMissingFeatures("ROTATION")
class Rotation(ut.TestCase):
    s = espressomd.System(box_l=[1.0, 1.0, 1.0])
    s.cell_system.skin = 0
    s.time_step = 0.01

    def tearDown(self):
        self.s.part.clear()

    def test_langevin(self):
        """Applies langevin thermostat and checks that correct axes get
           thermalized"""
        s = self.s
        s.thermostat.set_langevin(gamma=1, kT=1, seed=42)
        for rot_x, rot_y, rot_z in itertools.product((False, True), repeat=3):
            p = s.part.add(pos=(0, 0, 0), rotation=(rot_x, rot_y, rot_z),
                           quat=(1, 0, 0, 0), omega_body=(0, 0, 0),
                           torque_lab=(0, 0, 0))
            s.integrator.run(500)
            self.validate(p, rot_x, 0)
            self.validate(p, rot_y, 1)
            self.validate(p, rot_z, 2)
            s.part.clear()

    def validate(self, p, rotate, coord):
        if rotate:
            self.assertNotEqual(p.torque_lab[coord], 0)
            self.assertNotEqual(p.omega_body[coord], 0)
        else:
            # self.assertEqual(p.torque_lab[coord], 0)
            self.assertEqual(p.omega_body[coord], 0)

    @utx.skipIfMissingFeatures("EXTERNAL_FORCES")
    def test_axes_changes(self):
        """Verifies that rotation axes in body and space frame stay the same
           and other axes don't"""
        s = self.s
        p = s.part.add(pos=(0.9, 0.9, 0.9), ext_torque=(1, 1, 1))
        s.thermostat.turn_off()
        for direction in (0, 1, 2):
            # Reset orientation
            p.quat = [1, 0, 0, 0]

            # Enable rotation in a single direction
            rot = [0, 0, 0]
            rot[direction] = 1
            p.rotation = rot

            s.integrator.run(130)

            # Check other axes:
            for axis in [1, 0, 0], [0, 1, 0], [0, 0, 1]:
                if rot == axis:
                    # The axis for which rotation is on should coincide in body
                    # and space frame
                    self.assertAlmostEqual(
                        np.dot(rot, p.convert_vector_body_to_space(rot)), 1, places=8)
                else:
                    # For non-rotation axis, body and space frame should differ
                    self.assertLess(
                        np.dot(axis, p.convert_vector_body_to_space(axis)), 0.95)

    def test_frame_conversion_and_rotation(self):
        s = self.s
        p = s.part.add(pos=np.random.random(3), rotation=(1, 1, 1))

        # Space and body frame co-incide?
        np.testing.assert_allclose(
            np.copy(p.director), p.convert_vector_body_to_space((0, 0, 1)), atol=1E-10)

        # Random vector should still co-incide
        v = (1., 5.5, 17)
        np.testing.assert_allclose(
            v, p.convert_vector_space_to_body(v), atol=1E-10)
        np.testing.assert_allclose(
            v, p.convert_vector_body_to_space(v), atol=1E-10)

        # Particle rotation

        p.rotate((1, 2, 0), np.pi / 4)
        # Check angle for director
        self.assertAlmostEqual(
            np.arccos(np.dot(p.director, (0, 0, 1))), np.pi / 4, delta=1E-10)
        # Check other vector
        v = (5, -7, 3)
        v_r = p.convert_vector_body_to_space(v)
        self.assertAlmostEqual(np.dot(v, v), np.dot(v_r, v_r), delta=1e-10)
        np.testing.assert_allclose(
            p.convert_vector_space_to_body(v_r), v, atol=1E-10)

        # Rotation axis should co-incide
        np.testing.assert_allclose(
            (1, 2, 0), p.convert_vector_body_to_space((1, 2, 0)))

        # Check rotation axis with all elements set
        p.rotate(axis=(-5, 2, 17), angle=1.)
        v = (5, -7, 3)
        v_r = p.convert_vector_body_to_space(v)
        self.assertAlmostEqual(np.dot(v, v), np.dot(v_r, v_r), delta=1e-10)
        np.testing.assert_allclose(
            p.convert_vector_space_to_body(v_r), v, atol=1E-10)

    def test_rotation_mpi_communication(self):
        s = self.s
        s.part.clear()
        # place particle in cell with MPI rank 0
        p = s.part.add(pos=0.01 * self.s.box_l, rotation=(1, 1, 1))
        p.rotate((1, 0, 0), -np.pi / 2)
        np.testing.assert_array_almost_equal(
            np.copy(p.director), [0, 1, 0], decimal=10)
        # place particle in cell with MPI rank N-1
        p = s.part.add(pos=0.99 * self.s.box_l, rotation=(1, 1, 1))
        p.rotate((1, 0, 0), -np.pi / 2)
        np.testing.assert_array_almost_equal(
            np.copy(p.director), [0, 1, 0], decimal=10)


if __name__ == "__main__":
    ut.main()
