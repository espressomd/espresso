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
import sys
import unittest as ut
import unittest_decorators as utx
import numpy as np
import math

import espressomd
import espressomd.interactions
import espressomd.shapes
import tests_common


@utx.skipIfMissingFeatures(["LENNARD_JONES_GENERIC"])
class ShapeBasedConstraintTest(ut.TestCase):

    box_l = 30.
    system = espressomd.System(box_l=3 * [box_l])

    def tearDown(self):
        self.system.part.clear()
        self.system.constraints.clear()

    def pos_on_surface(self, theta, v, semiaxis0, semiaxis1,
                       semiaxis2, center=np.array([15, 15, 15])):
        """Return position on ellipsoid surface."""
        pos = np.array([semiaxis0 * np.sqrt(1. - v**2) * np.cos(theta),
                        semiaxis1 * np.sqrt(1. - v**2) * np.sin(theta),
                        semiaxis2 * v])
        return pos + center

    def test_hollow_conical_frustum(self):
        """
        Test implementation of conical frustum shape.

        """
        R1 = 5.0
        R2 = 10.0
        LENGTH = 15.0
        D = 2.4

        def z(y, r1, r2, l): return l / (r1 - r2) * \
            y + l / 2. - l * r1 / (r1 - r2)

        shape = espressomd.shapes.HollowConicalFrustum(center=[0.0, 0.0, 0.0], axis=[
            0, 0, 1], r1=R1, r2=R2, thickness=0.0, length=LENGTH)
        y_vals = np.linspace(R1, R2, 100)
        for y in y_vals:
            dist = shape.calc_distance(position=[0.0, y, z(y, R1, R2, LENGTH)])
            self.assertAlmostEqual(dist[0], 0.0)

        shape = espressomd.shapes.HollowConicalFrustum(center=[0.0, 0.0, 0.0], axis=[
            0, 0, 1], r1=R1, r2=R2, thickness=D, length=LENGTH, direction=-1)
        for y in y_vals:
            dist = shape.calc_distance(position=[0.0, y, z(y, R1, R2, LENGTH)])
            self.assertAlmostEqual(dist[0], 0.5 * D)

        np.testing.assert_almost_equal(np.copy(shape.center), [0.0, 0.0, 0.0])
        np.testing.assert_almost_equal(np.copy(shape.axis), [0, 0, 1])
        self.assertEqual(shape.r1, R1)
        self.assertEqual(shape.r2, R2)
        self.assertEqual(shape.thickness, D)
        self.assertEqual(shape.length, LENGTH)
        self.assertEqual(shape.direction, -1)

        shape = espressomd.shapes.HollowConicalFrustum(center=[0.0, 0.0, 0.0], axis=[
            0, 0, 1], r1=R1, r2=R2, thickness=D, length=LENGTH)
        for y in y_vals:
            dist = shape.calc_distance(position=[0.0, y, z(y, R1, R2, LENGTH)])
            self.assertAlmostEqual(dist[0], -0.5 * D)

        # check sign of dist
        shape = espressomd.shapes.HollowConicalFrustum(center=[0.0, 0.0, 0.0], axis=[
            0, 0, 1], r1=R1, r2=R1, thickness=D, length=LENGTH)
        self.assertLess(shape.calc_distance(
            position=[0.0, R1, 0.25 * LENGTH])[0], 0.0)
        self.assertLess(shape.calc_distance(
            position=[0.0, R1 + (0.5 - sys.float_info.epsilon) * D, 0.25 * LENGTH])[0], 0.0)
        self.assertGreater(shape.calc_distance(
            position=[0.0, R1 + (0.5 + sys.float_info.epsilon) * D, 0.25 * LENGTH])[0], 0.0)
        self.assertGreater(shape.calc_distance(
            position=[0.0, R1 - (0.5 + sys.float_info.epsilon) * D, 0.25 * LENGTH])[0], 0.0)

        shape = espressomd.shapes.HollowConicalFrustum(center=[0.0, 0.0, 0.0], axis=[
            0, 0, 1], r1=R1, r2=R1, thickness=D, length=LENGTH, direction=-1)
        self.assertGreater(shape.calc_distance(
            position=[0.0, R1, 0.25 * LENGTH])[0], 0.0)
        self.assertGreater(shape.calc_distance(
            position=[0.0, R1 + (0.5 - sys.float_info.epsilon) * D, 0.25 * LENGTH])[0], 0.0)
        self.assertLess(shape.calc_distance(
            position=[0.0, R1 + (0.5 + sys.float_info.epsilon) * D, 0.25 * LENGTH])[0], 0.0)
        self.assertLess(shape.calc_distance(
            position=[0.0, R1 - (0.5 + sys.float_info.epsilon) * D, 0.25 * LENGTH])[0], 0.0)

    def test_sphere(self):
        """Checks geometry of an inverted sphere

        """
        rad = self.box_l / 2.0
        sphere_shape = espressomd.shapes.Sphere(
            center=3 * [rad],
            radius=rad,
            direction=-1)
        phi_steps = 11
        theta_steps = 11
        for distance in {-1.2, 2.6}:
            for phi in range(phi_steps):
                phi_angle = phi / phi_steps * 2.0 * math.pi
                for theta in range(theta_steps):
                    theta_angle = theta / theta_steps * math.pi
                    pos = np.array(
                        [math.cos(phi_angle) * math.sin(theta_angle)
                         * (rad + distance),
                         math.sin(phi_angle) * math.sin(theta_angle)
                         * (rad + distance),
                         math.cos(theta_angle) * (rad + distance)]) + rad

                    shape_dist, _ = sphere_shape.calc_distance(
                        position=pos.tolist())
                    self.assertAlmostEqual(shape_dist, -distance)

    def test_ellipsoid(self):
        """Checks that distance of particles on the ellipsoid constraint's surface is zero.
        For the case of a spherical ellipsoid, also several non-zero distances are tested.

        """
        system = self.system
        system.time_step = 0.01
        system.cell_system.skin = 0.4
        system.part.add(pos=[0., 0., 0.], type=0)

        # abuse generic LJ to measure distance via the potential V(r) = r
        system.non_bonded_inter[0, 1].generic_lennard_jones.set_params(
            epsilon=1., sigma=1., cutoff=7., shift=0., offset=0., e1=-1, e2=0, b1=1., b2=0.)

        N = 10

        # check oblate ellipsoid

        semiaxes = [2.18, 5.45]
        e = espressomd.shapes.Ellipsoid(
            a=semiaxes[0],
            b=semiaxes[1],
            center=3 * [self.box_l / 2.],
            direction=+1)

        constraint_e = espressomd.constraints.ShapeBasedConstraint(
            shape=e, particle_type=1, penetrable=True)
        const1 = system.constraints.add(constraint_e)

        for i in range(N):
            for j in range(N):
                theta = 2. * i / float(N) * np.pi
                v = j / float(N - 1) * 2. - 1
                pos = self.pos_on_surface(
                    theta, v, semiaxes[0], semiaxes[1], semiaxes[1])
                system.part[0].pos = pos
                system.integrator.run(recalc_forces=True, steps=0)
                energy = system.analysis.energy()
                self.assertAlmostEqual(energy["total"], 0., places=6)

        system.constraints.remove(const1)

        # check prolate ellipsoid

        semiaxes = [3.61, 2.23]
        e = espressomd.shapes.Ellipsoid(
            a=semiaxes[0],
            b=semiaxes[1],
            center=3 * [self.box_l / 2.],
            direction=+1)

        constraint_e = espressomd.constraints.ShapeBasedConstraint(
            shape=e, particle_type=1, penetrable=True)
        const1 = system.constraints.add(constraint_e)

        for i in range(N):
            for j in range(N):
                theta = 2. * i / float(N) * np.pi
                v = j / float(N - 1) * 2. - 1
                pos = self.pos_on_surface(
                    theta, v, semiaxes[0], semiaxes[1], semiaxes[1])
                system.part[0].pos = pos
                system.integrator.run(recalc_forces=True, steps=0)
                energy = system.analysis.energy()
                self.assertAlmostEqual(energy["total"], 0., places=6)

        # check sphere (multiple distances from surface)

        # change ellipsoid parameters instead of creating a new constraint
        e.a = 1.
        e.b = 1.

        radii = np.linspace(1., 6.5, 7)

        for i in range(N):
            for j in range(N):
                theta = 2. * i / float(N) * np.pi
                v = j / float(N - 1) * 2. - 1
                for r in radii:
                    pos = self.pos_on_surface(theta, v, r, r, r)
                    system.part[0].pos = pos
                    system.integrator.run(recalc_forces=True, steps=0)
                    energy = system.analysis.energy()
                    self.assertAlmostEqual(energy["total"], r - 1.)
        # Reset the interaction to zero
        system.non_bonded_inter[0, 1].generic_lennard_jones.set_params(
            epsilon=0., sigma=0., cutoff=0., shift=0., offset=0., e1=0, e2=0, b1=0., b2=0.)

    def test_cylinder(self):
        """Tests if shape based constraints can be added to a system both by
        (1) defining a constraint object which is then added
        (2) and via keyword arguments.
        Checks that cylinder constraints with LJ interactions exert forces
        on a test particle (that is, the constraints do what they should).

        """
        system = self.system
        system.time_step = 0.01
        system.cell_system.skin = 0.4
        rad = self.box_l / 2.0
        length = self.box_l / 2.0

        system.part.add(id=0, pos=[rad, 1.02, rad], type=0)

        # check force calculation of a cylinder without top and bottom
        interaction_dir = -1  # constraint is directed inwards
        cylinder_shape = espressomd.shapes.Cylinder(
            center=3 * [rad],
            axis=[0, 0, 1],
            direction=interaction_dir,
            radius=rad,
            length=self.box_l + 5)  # +5 in order to have no top or bottom
        penetrability = False  # impenetrable
        outer_cylinder_constraint = espressomd.constraints.ShapeBasedConstraint(
            shape=cylinder_shape, particle_type=1, penetrable=penetrability)
        outer_cylinder_wall = system.constraints.add(outer_cylinder_constraint)
        system.non_bonded_inter[0, 1].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0, cutoff=2.0, shift=0)
        system.integrator.run(0)  # update forces

        self.assertAlmostEqual(outer_cylinder_constraint.min_dist(), 1.02)

        # test summed forces on cylinder wall
        self.assertAlmostEqual(
            -1.0 * outer_cylinder_wall.total_force()[1],
            tests_common.lj_force(
                espressomd,
                cutoff=2.0,
                offset=0.,
                eps=1.0,
                sig=1.0,
                r=1.02),
            places=10)  # minus for Newton's third law

        # check whether total_summed_outer_normal_force is correct
        y_part2 = self.box_l - 1.02
        system.part.add(id=1, pos=[rad, y_part2, rad], type=0)
        system.integrator.run(0)

        dist_part2 = self.box_l - y_part2
        self.assertAlmostEqual(outer_cylinder_wall.total_force()[2], 0.0)
        self.assertAlmostEqual(
            outer_cylinder_wall.total_normal_force(),
            2 *
            tests_common.lj_force(
                espressomd,
                cutoff=2.0,
                offset=0.,
                eps=1.0,
                sig=1.0,
                r=dist_part2))

        # Test the geometry of a cylinder with top and bottom
        cylinder_shape_finite = espressomd.shapes.Cylinder(
            center=3 * [rad],
            axis=[0, 0, 1],
            direction=1,
            radius=rad,
            length=length)

        phi_steps = 11
        for distance in {-3.6, 2.8}:
            for z in range(int(self.box_l)):
                center = np.array([rad, rad, z])
                start_point = np.array([rad, 2 * rad - distance, z])
                for phi in range(phi_steps):
                    # Rotation around the axis of the cylinder
                    phi_angle = phi / phi_steps * 2.0 * math.pi
                    phi_rot_matrix = np.array(
                        [[math.cos(phi_angle), -math.sin(phi_angle), 0.0],
                         [math.sin(phi_angle), math.cos(phi_angle), 0.0],
                         [0.0, 0.0, 1.0]])
                    phi_rot_point = np.dot(
                        phi_rot_matrix, start_point - center) + center

                    shape_dist, _ = cylinder_shape_finite.calc_distance(
                        position=phi_rot_point.tolist())

                    dist = -distance
                    if distance > 0.0:
                        if z < (self.box_l - length) / 2.0 + distance:
                            dist = (self.box_l - length) / 2.0 - z
                        elif z > (self.box_l + length) / 2.0 - distance:
                            dist = z - (self.box_l + length) / 2.0
                        else:
                            dist = -distance
                    else:
                        if z < (self.box_l - length) / 2.0:
                            z_dist = (self.box_l - length) / 2.0 - z
                            dist = math.sqrt(z_dist**2 + distance**2)
                        elif z > (self.box_l + length) / 2.0:
                            z_dist = z - (self.box_l + length) / 2.0
                            dist = math.sqrt(z_dist**2 + distance**2)
                        else:
                            dist = -distance

                    self.assertAlmostEqual(shape_dist, dist)
        # Reset
        system.non_bonded_inter[0, 1].lennard_jones.set_params(
            epsilon=0.0, sigma=0.0, cutoff=0.0, shift=0)

    def test_spherocylinder(self):
        """Checks that spherocylinder constraints with LJ interactions exert
        forces on a test particle (that is, the constraints do what they should)
        using geometrical parameters of (1) an infinite cylinder and (2) a
        finite spherocylinder.

        """
        system = self.system
        system.time_step = 0.01
        system.cell_system.skin = 0.4

        system.part.add(
            id=0, pos=[self.box_l / 2.0, 1.02, self.box_l / 2.0], type=0)

        # check force calculation of spherocylinder constraint
        # (1) infinite cylinder
        interaction_dir = -1  # constraint is directed inwards
        spherocylinder_shape = espressomd.shapes.SpheroCylinder(
            center=3 * [self.box_l / 2.0],
            axis=[0, 0, 1],
            direction=interaction_dir,
            radius=self.box_l / 2.0,
            length=self.box_l + 5)  # +5 in order to have no top or bottom
        penetrability = False  # impenetrable
        outer_cylinder_constraint = espressomd.constraints.ShapeBasedConstraint(
            shape=spherocylinder_shape, particle_type=1, penetrable=penetrability)
        system.constraints.add(outer_cylinder_constraint)
        system.non_bonded_inter[0, 1].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0, cutoff=2.0, shift=0)
        system.integrator.run(0)  # update forces

        self.assertAlmostEqual(outer_cylinder_constraint.min_dist(), 1.02)

        # test summed forces on cylinder wall
        self.assertAlmostEqual(
            -1.0 * outer_cylinder_constraint.total_force()[1],
            tests_common.lj_force(
                espressomd,
                cutoff=2.0,
                offset=0.,
                eps=1.0,
                sig=1.0,
                r=1.02),
            places=10)  # minus for Newton's third law

        # check whether total_summed_outer_normal_force is correct
        y_part2 = self.box_l - 1.02
        system.part.add(
            id=1, pos=[self.box_l / 2.0, y_part2, self.box_l / 2.0], type=0)
        system.integrator.run(0)

        dist_part2 = self.box_l - y_part2
        self.assertAlmostEqual(outer_cylinder_constraint.total_force()[2], 0.0)
        self.assertAlmostEqual(outer_cylinder_constraint.total_normal_force(),
                               2 * tests_common.lj_force(
                                   espressomd, cutoff=2.0, offset=0.,
                                   eps=1.0, sig=1.0, r=dist_part2))

        # Reset
        system.part.clear()
        system.constraints.clear()
        system.non_bonded_inter[0, 1].lennard_jones.set_params(
            epsilon=0.0, sigma=0.0, cutoff=0.0, shift=0)

        # (2) finite spherocylinder
        system.part.clear()
        interaction_dir = -1  # constraint is directed inwards
        spherocylinder_shape = espressomd.shapes.SpheroCylinder(
            center=3 * [self.box_l / 2.0],
            axis=[0, 1, 0],
            direction=interaction_dir,
            radius=10.0,
            length=6.0)
        penetrability = True  # penetrable
        inner_cylinder_constraint = espressomd.constraints.ShapeBasedConstraint(
            shape=spherocylinder_shape, particle_type=1, penetrable=penetrability)
        system.constraints.add(inner_cylinder_constraint)

        # V(r) = r
        system.non_bonded_inter[0, 1].generic_lennard_jones.set_params(
            epsilon=1., sigma=1., cutoff=10., shift=0., offset=0., e1=-1, e2=0, b1=1., b2=0.)

        # check hemispherical caps (multiple distances from surface)
        N = 10
        radii = np.linspace(1., 10., 10)
        system.part.add(pos=[0., 0., 0.], type=0)
        for i in range(6):
            for j in range(N):
                theta = 2. * i / float(N) * np.pi
                v = j / float(N - 1) * 2. - 1
                for r in radii:
                    pos = self.pos_on_surface(theta, v, r, r, r) + [0, 3, 0]
                    system.part[0].pos = pos
                    system.integrator.run(recalc_forces=True, steps=0)
                    energy = system.analysis.energy()
                    self.assertAlmostEqual(energy["total"], 10. - r)

        # Reset
        system.non_bonded_inter[0, 1].generic_lennard_jones.set_params(
            epsilon=0., sigma=0., cutoff=0., shift=0., offset=0., e1=0, e2=0, b1=0., b2=0.)

    def test_wall_forces(self):
        """Tests if shape based constraints can be added to a system both by
        (1) defining a constraint object which is then added
        (2) and via keyword arguments.
        Checks that wall constraints with LJ interactions exert forces
        on a test particle (that is, the constraints do what they should).

        """
        system = self.system
        system.time_step = 0.01
        system.part.add(id=0, pos=[5., 1.21, 0.83], type=0)

        # Check forces are initialized to zero
        f_part = system.part[0].f

        self.assertEqual(f_part[0], 0.)
        self.assertEqual(f_part[1], 0.)
        self.assertEqual(f_part[2], 0.)

        system.non_bonded_inter[0, 1].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0, cutoff=2.0, shift=0)
        system.non_bonded_inter[0, 2].lennard_jones.set_params(
            epsilon=1.5, sigma=1.0, cutoff=2.0, shift=0)

        shape_xz = espressomd.shapes.Wall(normal=[0., 1., 0.], dist=0.)
        shape_xy = espressomd.shapes.Wall(normal=[0., 0., 1.], dist=0.)

        # (1)
        constraint_xz = espressomd.constraints.ShapeBasedConstraint(
            shape=shape_xz, particle_type=1)
        wall_xz = system.constraints.add(constraint_xz)

        # (2)
        wall_xy = system.constraints.add(shape=shape_xy, particle_type=2)

        system.integrator.run(0)  # update forces
        f_part = system.part[0].f

        self.assertEqual(f_part[0], 0.)
        self.assertAlmostEqual(
            f_part[1],
            tests_common.lj_force(
                espressomd,
                cutoff=2.0,
                offset=0.,
                eps=1.0,
                sig=1.0,
                r=1.21),
            places=10)
        self.assertAlmostEqual(
            f_part[2],
            tests_common.lj_force(
                espressomd,
                cutoff=2.0,
                offset=0.,
                eps=1.5,
                sig=1.0,
                r=0.83),
            places=10)

        # test summed forces on walls
        self.assertAlmostEqual(
            -1.0 * wall_xz.total_force()[1],
            tests_common.lj_force(
                espressomd,
                cutoff=2.0,
                offset=0.,
                eps=1.0,
                sig=1.0,
                r=1.21),
            places=10)  # minus for Newton's third law
        self.assertAlmostEqual(
            -1.0 * wall_xy.total_force()[2],
            tests_common.lj_force(
                espressomd,
                cutoff=2.0,
                offset=0.,
                eps=1.5,
                sig=1.0,
                r=0.83),
            places=10)

        # check whether total_normal_force is correct
        self.assertAlmostEqual(
            wall_xy.total_normal_force(),
            tests_common.lj_force(
                espressomd,
                cutoff=2.0,
                offset=0.,
                eps=1.5,
                sig=1.0,
                r=0.83),
            places=10)

        # this one is closer and should get the mindist()
        system.part.add(pos=[5., 1.20, 0.82], type=0)
        self.assertAlmostEqual(constraint_xz.min_dist(), system.part[1].pos[1])
        self.assertAlmostEqual(wall_xz.min_dist(), system.part[1].pos[1])
        self.assertAlmostEqual(wall_xy.min_dist(), system.part[1].pos[2])

        # Reset
        system.non_bonded_inter[0, 1].lennard_jones.set_params(
            epsilon=0.0, sigma=0.0, cutoff=0.0, shift=0)
        system.non_bonded_inter[0, 2].lennard_jones.set_params(
            epsilon=0.0, sigma=0.0, cutoff=0.0, shift=0)

    def test_slitpore(self):
        """Checks that slitpore constraints with LJ interactions exert forces
        on a test particle (that is, the constraints do what they should).

        """
        system = self.system
        system.time_step = 0.01
        system.cell_system.skin = 0.4

        # check force calculation of slitpore constraint
        slitpore_shape = espressomd.shapes.Slitpore(
            channel_width=5,
            lower_smoothing_radius=2,
            upper_smoothing_radius=3,
            pore_length=15,
            pore_mouth=20,
            pore_width=10,
            dividing_plane=self.box_l / 2)
        slitpore_constraint = espressomd.constraints.ShapeBasedConstraint(
            shape=slitpore_shape, particle_type=1, penetrable=True)
        system.constraints.add(slitpore_constraint)
        # V(r) = r
        system.non_bonded_inter[0, 1].generic_lennard_jones.set_params(
            epsilon=1., sigma=1., cutoff=10., shift=0., offset=0., e1=-1, e2=0, b1=1., b2=0.)

        system.part.add(pos=[0., 0., 0.], type=0)
        x = self.box_l / 2.0
        parameters = [
            ([x, x, 1.], -4., [0., 0., -1.]),  # outside channel
            ([x, x, 15.], 5., [-1., 0., 0.]),  # inside channel
            ([x, x, 5.], 0., [0., 0., 0.]),  # on channel bottom surface
            ([x + 5., x, 15.], 0., [0., 0., 0.]),  # on channel side surface
            ([x, x, 15.], 5., [-1., 0., 0.]),  # within mouth
            ([x, x, 25.], 0., [0., 0., 0.]),  # on wall surface
            ([x, x, 27.], -2., [0., 0., 1.]),  # outside wall
        ]
        for pos, ref_mindist, ref_force in parameters:
            system.part[0].pos = pos
            system.integrator.run(recalc_forces=True, steps=0)
            self.assertEqual(slitpore_constraint.min_dist(), ref_mindist)
            np.testing.assert_almost_equal(
                np.copy(slitpore_constraint.total_force()), ref_force, 10)

        # Reset
        system.non_bonded_inter[0, 1].generic_lennard_jones.set_params(
            epsilon=0., sigma=0., cutoff=0., shift=0., offset=0., e1=0, e2=0, b1=0., b2=0.)

    def test_rhomboid(self):
        """Checks that rhomboid constraints with LJ interactions exert forces
        on a test particle (that is, the constraints do what they should)
        using the geometrical parameters of (1) a cuboid and (2) a rhomboid.

        """
        system = self.system
        system.time_step = 0.01
        system.cell_system.skin = 0.4

        # check force calculation of rhomboid constraint
        # (1) using a cuboid
        interaction_dir = +1  # constraint is directed outwards
        length = np.array([-5.0, 6.0, 7.0])  # dimension of the cuboid
        corner = np.array(3 * [self.box_l / 2.0])
        rhomboid_shape = espressomd.shapes.Rhomboid(
            corner=corner,
            a=[length[0], 0.0, 0.0],  # cube
            b=[0.0, length[1], 0.0],
            c=[0.0, 0.0, length[2]],
            direction=interaction_dir
        )
        penetrability = False  # impenetrable
        rhomboid_constraint = espressomd.constraints.ShapeBasedConstraint(
            shape=rhomboid_shape, particle_type=1, penetrable=penetrability)
        rhomboid_constraint = system.constraints.add(rhomboid_constraint)

        system.non_bonded_inter[0, 1].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0, cutoff=2.0, shift=0)
        system.part.add(id=0, pos=[self.box_l / 2.0 + length[0] / 2.0,
                                   self.box_l / 2.0 + length[1] / 2.0,
                                   self.box_l / 2.0 - 1], type=0)
        system.integrator.run(0)  # update forces
        f_part = system.part[0].f
        self.assertEqual(rhomboid_constraint.min_dist(), 1.)
        self.assertEqual(f_part[0], 0.)
        self.assertEqual(f_part[1], 0.)
        self.assertAlmostEqual(
            -f_part[2],
            tests_common.lj_force(
                espressomd,
                cutoff=2.,
                offset=0.,
                eps=1.,
                sig=1.,
                r=1.),
            places=10)
        self.assertAlmostEqual(
            rhomboid_constraint.total_normal_force(),
            tests_common.lj_force(
                espressomd,
                cutoff=2.,
                offset=0.,
                eps=1.,
                sig=1.,
                r=1.),
            places=10)

        x_range = 12
        y_range = 12
        z_range = 12
        for x in range(x_range):
            for y in range(y_range):
                for z in range(z_range):
                    pos = np.array(
                        [x + (self.box_l + length[0] - x_range) / 2.0,
                         y + (self.box_l + length[1] - y_range) / 2.0,
                         z + (self.box_l + length[2] - z_range) / 2.0])
                    shape_dist, shape_dist_vec = rhomboid_shape.calc_distance(
                        position=pos.tolist())

                    outside = False
                    edge_case = False
                    dist_vec = np.array([0.0, 0.0, 0.0])

                    # check if outside or inside
                    if(pos[0] < (self.box_l + length[0] - abs(length[0])) / 2.0 or
                       pos[0] > (self.box_l + length[0] + abs(length[0])) / 2.0 or
                       pos[1] < (self.box_l + length[1] - abs(length[1])) / 2.0 or
                       pos[1] > (self.box_l + length[1] + abs(length[1])) / 2.0 or
                       pos[2] < (self.box_l + length[2] - abs(length[2])) / 2.0 or
                       pos[2] > (self.box_l + length[2] + abs(length[2])) / 2.0):
                        outside = True

                    if outside:
                        for i in range(3):
                            if pos[i] < (self.box_l + length[i] -
                                         abs(length[i])) / 2.0:
                                dist_vec[i] = pos[i] - (
                                    self.box_l + length[i] - abs(length[i])) / 2.0
                            elif pos[i] > (self.box_l + length[i] + abs(length[i])) / 2.0:
                                dist_vec[i] = pos[i] - (
                                    self.box_l + length[i] + abs(length[i])) / 2.0
                            else:
                                dist_vec[i] = 0.0
                        dist = np.linalg.norm(dist_vec)
                    else:
                        dist = self.box_l
                        c1 = pos - corner
                        c2 = corner + length - pos
                        abs_c1c2 = np.abs(np.concatenate((c1, c2)))
                        dist = np.amin(abs_c1c2)
                        where = np.argwhere(dist == abs_c1c2)
                        if len(where) > 1:
                            edge_case = True
                        for which in where:
                            if which < 3:
                                dist_vec[which] = dist * np.sign(c1[which])
                            else:
                                dist_vec[which - 3] = -dist * \
                                    np.sign(c2[which - 3])

                        dist *= -interaction_dir

                    if edge_case:
                        for i in range(3):
                            if shape_dist_vec[i] != 0.0:
                                self.assertAlmostEqual(
                                    abs(shape_dist_vec[i]), abs(dist_vec[i]))
                    else:
                        self.assertAlmostEqual(shape_dist_vec[0], dist_vec[0])
                        self.assertAlmostEqual(shape_dist_vec[1], dist_vec[1])
                        self.assertAlmostEqual(shape_dist_vec[2], dist_vec[2])
                    self.assertAlmostEqual(shape_dist, dist)

        # (2) using a rhomboid
        rhomboid_shape.a = [5., 5., 0.]  # rhomboid
        rhomboid_shape.b = [0., 0., 5.]
        rhomboid_shape.c = [0., 5., 0.]

        system.part[0].pos = [self.box_l / 2.0 + 2.5,
                              self.box_l / 2.0 + 2.5,
                              self.box_l / 2.0 - 1]

        system.integrator.run(0)  # update forces
        self.assertEqual(rhomboid_constraint.min_dist(), 1.)
        self.assertAlmostEqual(
            rhomboid_constraint.total_normal_force(),
            tests_common.lj_force(
                espressomd,
                cutoff=2.,
                offset=0.,
                eps=1.,
                sig=1.,
                r=1.),
            places=10)

        system.part[0].pos = system.part[0].pos - [0., 1., 0.]
        system.integrator.run(0)  # update forces
        self.assertAlmostEqual(
            rhomboid_constraint.min_dist(), 1.2247448714, 10)
        self.assertAlmostEqual(
            rhomboid_constraint.total_normal_force(),
            tests_common.lj_force(
                espressomd,
                cutoff=2.,
                offset=0.,
                eps=1.,
                sig=1.,
                r=1.2247448714),
            places=10)

        # Reset
        system.non_bonded_inter[0, 1].lennard_jones.set_params(
            epsilon=0.0, sigma=0.0, cutoff=0.0, shift=0)

    def test_torus(self):
        """Checks that torus constraints with LJ interactions exert forces
        on a test particle (that is, the constraints do what they should).

        """
        system = self.system
        system.time_step = 0.01
        system.cell_system.skin = 0.4

        interaction_dir = 1  # constraint is directed inwards
        radius = self.box_l / 4.0
        tube_radius = self.box_l / 6.0
        part_offset = 1.2

        system.part.add(
            id=0, pos=[self.box_l / 2.0, self.box_l / 2.0 + part_offset, self.box_l / 2.0], type=0)

        # check force calculation of cylinder constraint
        torus_shape = espressomd.shapes.Torus(
            center=3 * [self.box_l / 2.0],
            normal=[0, 0, 1],
            direction=interaction_dir,
            radius=radius,
            tube_radius=tube_radius)
        penetrability = False  # impenetrable
        torus_constraint = espressomd.constraints.ShapeBasedConstraint(
            shape=torus_shape, particle_type=1, penetrable=penetrability)
        torus_wall = system.constraints.add(torus_constraint)
        system.non_bonded_inter[0, 1].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0, cutoff=2.0, shift=0)
        system.integrator.run(0)  # update forces

        self.assertAlmostEqual(torus_constraint.min_dist(),
                               radius - tube_radius - part_offset)

        # test summed forces on torus wall
        self.assertAlmostEqual(
            -1.0 * torus_wall.total_force()[1],
            tests_common.lj_force(
                espressomd,
                cutoff=2.0,
                offset=0.,
                eps=1.0,
                sig=1.0,
                r=torus_constraint.min_dist()),
            places=10)  # minus for Newton's third law

        # check whether total_summed_outer_normal_force is correct
        y_part2 = self.box_l / 2.0 + 2.0 * radius - part_offset
        system.part.add(
            id=1, pos=[self.box_l / 2.0, y_part2, self.box_l / 2.0], type=0)
        system.integrator.run(0)

        self.assertAlmostEqual(torus_wall.total_force()[1], 0.0)
        self.assertAlmostEqual(torus_wall.total_normal_force(), 2 * tests_common.lj_force(
            espressomd, cutoff=2.0, offset=0., eps=1.0, sig=1.0,
            r=radius - tube_radius - part_offset))

        # Test the geometry of the shape directly
        phi_steps = 11
        theta_steps = 11
        center = np.array([self.box_l / 2.0,
                           self.box_l / 2.0,
                           self.box_l / 2.0])
        tube_center = np.array([self.box_l / 2.0,
                                self.box_l / 2.0 + radius,
                                self.box_l / 2.0])

        for distance in {1.02, -0.7}:
            start_point = np.array([self.box_l / 2.0,
                                    self.box_l / 2.0 + radius -
                                    tube_radius - distance,
                                    self.box_l / 2.0])
            for phi in range(phi_steps):
                for theta in range(theta_steps):
                    # Rotation around the tube
                    theta_angle = theta / theta_steps * 2.0 * math.pi
                    theta_rot_matrix = np.array(
                        [[1.0, 0.0, 0.0],
                         [0.0, math.cos(theta_angle), -math.sin(theta_angle)],
                         [0.0, math.sin(theta_angle), math.cos(theta_angle)]])
                    theta_rot_point = np.dot(
                        theta_rot_matrix,
                        start_point - tube_center)
                    theta_rot_point += tube_center

                    # Rotation around the center of the torus
                    phi_angle = phi / phi_steps * 2.0 * math.pi
                    phi_rot_matrix = np.array(
                        [[math.cos(phi_angle), -math.sin(phi_angle), 0.0],
                         [math.sin(phi_angle), math.cos(phi_angle), 0.0],
                         [0.0, 0.0, 1.0]])
                    phi_rot_point = np.dot(
                        phi_rot_matrix,
                        theta_rot_point - center) + center

                    shape_dist, _ = torus_shape.calc_distance(
                        position=phi_rot_point.tolist())
                    self.assertAlmostEqual(shape_dist, distance)

        # Reset
        system.non_bonded_inter[0, 1].lennard_jones.set_params(
            epsilon=0.0, sigma=0.0, cutoff=0.0, shift=0)


if __name__ == "__main__":
    ut.main()
