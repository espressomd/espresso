from __future__ import division, print_function

import unittest as ut

import espressomd

from espressomd import interactions
import espressomd.shapes as shapes
import tests_common


@ut.skipIf(not espressomd.has_features(["CONSTRAINTS", "LENNARD_JONES"]),
           "Features not available, skipping test!")
class ShapeBasedConstraintTest(ut.TestCase):
    """Tests if shape based constraints can be added to a system both by
    (1) defining a constraint object which is then added
    (2) and via keyword arguments.
    Checks that wall constraints with LJ interactions exert forces
    on a test particle (that is, the constraints do what they should).

    """

    box_l = 10.

    def prepare(self, system):
        system.box_l = [self.box_l, self.box_l, self.box_l]
        system.time_step = 0.01
        system.cell_system.skin = 0.4

        system.part.add(pos=[5., 1.21, 0.83], type=0)

        system.non_bonded_inter[0, 1].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0, cutoff=2.0, shift=0)
        system.non_bonded_inter[0, 2].lennard_jones.set_params(
            epsilon=1.5, sigma=1.0, cutoff=2.0, shift=0)


    def test(self):
        system = espressomd.System(box_l=[1.0, 1.0, 1.0])
        self.prepare(system)

        wy = shapes.Wall(normal=[0., 1., 0.], dist=0.)
        wz = shapes.Wall(normal=[0., 0., 1.], dist=0.)

        # (1)
        constraint_wy = espressomd.constraints.ShapeBasedConstraint(
            shape=wy, particle_type=1)
        wall_xy = system.constraints.add(constraint_wy)

        # (3)
        wall_xz = system.constraints.add(shape=wz, particle_type=2)

        # Check forces
        f_part = system.part[0].f


        self.assertEqual(f_part[0], 0.)
        self.assertEqual(f_part[1], 0.)
        self.assertEqual(f_part[2], 0.)

        system.integrator.run(0)  # update forces
        f_part = system.part[0].f

        self.assertEqual(f_part[0], 0.)
        self.assertAlmostEqual(f_part[1], tests_common.lj_force(espressomd, cutoff=2.0, offset=0.,
                                                                eps=1.0, sig=1.0, r=1.21), places=10)
        self.assertAlmostEqual(f_part[2], tests_common.lj_force(espressomd, cutoff=2.0, offset=0.,
                                                                eps=1.5, sig=1.0, r=0.83), places=10)

        # test forces on walls
        self.assertAlmostEqual(-1.0 * wall_xy.total_force()[1], tests_common.lj_force(espressomd, cutoff=2.0, offset=0.,
                                                                                      eps=1.0, sig=1.0, r=1.21), places=10)  # minus for newtons thrid law
        self.assertAlmostEqual(-1.0 * wall_xz.total_force()[2], tests_common.lj_force(espressomd, cutoff=2.0, offset=0.,
                                                                                      eps=1.5, sig=1.0, r=0.83), places=10)
        # this one is closer and should get the mindist()
        system.part.add(pos=[5., 1.20, 0.82], type=0)
        self.assertAlmostEqual(constraint_wy.min_dist(), system.part[1].pos[1])
        self.assertAlmostEqual(wall_xy.min_dist(), system.part[1].pos[1])
        self.assertAlmostEqual(wall_xz.min_dist(), system.part[1].pos[2])

        # Check removal
        for c in system.constraints:
            system.constraints.remove(c)

        for c in system.constraints:
            self.assertTrue(False)

if __name__ == "__main__":
    ut.main()
