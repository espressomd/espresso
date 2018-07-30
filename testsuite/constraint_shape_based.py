from __future__ import division, print_function

import unittest as ut
import numpy

import espressomd
import espressomd.interactions
import espressomd.shapes
import tests_common


@ut.skipIf(not espressomd.has_features(["LENNARD_JONES_GENERIC"]),
           "Features not available, skipping test!")
class ShapeBasedConstraintTest(ut.TestCase):

    box_l = 30.
    system = espressomd.System(box_l=3 * [box_l])

    def tearDown(self):
        self.system.part.clear()
        self.system.constraints.clear()

    def pos_on_surface(self, theta, v, semiaxis0, semiaxis1,
                       semiaxis2, center=numpy.array([15, 15, 15])):
        """Return postion on ellipsoid surface."""
        pos = numpy.array([semiaxis0 *
                           numpy.sqrt(1. -
                                      v *
                                      v) *
                           numpy.cos(theta), semiaxis1 *
                           numpy.sqrt(1. -
                                      v *
                                      v) *
                           numpy.sin(theta), semiaxis2 *
                           v])
        return pos + center

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
            center=[
                self.box_l / 2.,
                self.box_l / 2.,
                self.box_l / 2.],
            direction=+1)

        constraint_e = espressomd.constraints.ShapeBasedConstraint(
            shape=e, particle_type=1, penetrable=True)
        const1 = system.constraints.add(constraint_e)

        for i in range(N):
            for j in range(N):
                theta = 2. * i / float(N) * numpy.pi
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
            center=[
                self.box_l / 2.,
                self.box_l / 2.,
                self.box_l / 2.],
            direction=+1)

        constraint_e = espressomd.constraints.ShapeBasedConstraint(
            shape=e, particle_type=1, penetrable=True)
        const1 = system.constraints.add(constraint_e)

        for i in range(N):
            for j in range(N):
                theta = 2. * i / float(N) * numpy.pi
                v = j / float(N - 1) * 2. - 1
                pos = self.pos_on_surface(
                    theta, v, semiaxes[0], semiaxes[1], semiaxes[1])
                system.part[0].pos = pos
                system.integrator.run(recalc_forces=True, steps=0)
                energy = system.analysis.energy()
                self.assertAlmostEqual(energy["total"], 0., places=6)

        # check sphere (multiple distances from surface)

        # change ellipsoid paramters instead of creating a new constraint
        e.a = 1.
        e.b = 1.

        radii = numpy.linspace(1., 6.5, 7)

        for i in range(N):
            for j in range(N):
                theta = 2. * i / float(N) * numpy.pi
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

        system.part.add(id=0, pos=[self.box_l / 2.0, 1.02, self.box_l / 2.0], type=0)

        # check force calculation of cylinder constraint
        interactio_dir = -1  # constraint is directed inwards
        cylinder_shape = espressomd.shapes.Cylinder(
            center=[
                self.box_l /
                2.0,
                self.box_l /
                2.0,
                self.box_l /
                2.0],
            axis=[
                0,
                0,
                1],
            direction=interactio_dir,
            radius=self.box_l /
            2.0,
            length=self.box_l +
            5)  # +5 in order to have no top or bottom
        penetrability = False  # inpenetrable
        outer_cylinder_constraint = espressomd.constraints.ShapeBasedConstraint(
            shape=cylinder_shape, particle_type=1, penetrable=penetrability)
        outer_cylinder_wall = system.constraints.add(outer_cylinder_constraint)
        system.non_bonded_inter[0, 1].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0, cutoff=2.0, shift=0)
        system.non_bonded_inter[0, 2].lennard_jones.set_params(
            epsilon=1.5, sigma=1.0, cutoff=2.0, shift=0)
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
            places=10)  # minus for newtons thrid law
            
        #check whether total_summed_outer_normal_force is correct
        y_part2=self.box_l-1.02
        system.part.add(id=1, pos=[self.box_l / 2.0, y_part2, self.box_l / 2.0], type=0)
        system.integrator.run(0)
        
        dist_part2=self.box_l-y_part2
        self.assertAlmostEqual(outer_cylinder_wall.total_force()[2],0.0)
        self.assertAlmostEqual(outer_cylinder_wall.total_normal_force(),2*tests_common.lj_force( espressomd, cutoff=2.0, offset=0., eps=1.0, sig=1.0, r=dist_part2))
        
        # Reset
        system.non_bonded_inter[0, 1].lennard_jones.set_params(
            epsilon=0.0, sigma=0.0, cutoff=0.0, shift=0)
        system.non_bonded_inter[0, 2].lennard_jones.set_params(
            epsilon=0.0, sigma=0.0, cutoff=0.0, shift=0)

    def test_wall_forces(self):
        """Tests if shape based constraints can be added to a system both by
        (1) defining a constraint object which is then added
        (2) and via keyword arguments.
        Checks that wall constraints with LJ interactions exert forces
        on a test particle (that is, the constraints do what they should).

        """
        system = self.system
        system.part.add(id=0,pos=[5., 1.21, 0.83], type=0)

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

        # (3)
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
            places=10)  # minus for newtons thrid law
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

        #check whether total_normal_force is correct
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


if __name__ == "__main__":
    ut.main()
