# Tests if shape based constraints can be added to a system both by
#  (1) defining a constraint object which is then added
#  (2) and via keyword arguments.
# Checks, that (Wall-)constraints with LJ interactions exert forces
# on a test particle (that is, the constraints do what they should).

from __future__ import division, print_function

import unittest as ut

import espressomd

from espressomd import interactions
from espressomd.shapes import Wall


@ut.skipIf(not espressomd.has_features(["CONSTRAINTS", "LENNARD_JONES"]),
           "Features not available, skipping test!")
class ShapeBasedConstraintTest(ut.TestCase):

    def prepare(self, S):
        S.box_l = [10., 10., 10.]
        S.time_step = 0.01
        S.cell_system.skin = 0.4

        S.part.add(pos=[5., 1.21, 0.83], type=0)

        S.non_bonded_inter[0, 1].lennard_jones.set_params(
            epsilon=1.0, sigma=1.0, cutoff=2.0, shift=0)
        S.non_bonded_inter[0, 2].lennard_jones.set_params(
            epsilon=1.5, sigma=1.0, cutoff=2.0, shift=0)

    def lj_force(self, eps, sig, r):
        f_lj = 24.0 * eps * (2.0 * sig**12 / r**13 - sig**6 / r**7)
        return f_lj

    def test(self):
        S = espressomd.System()
        self.prepare(S)

        wy = Wall(normal=[0., 1., 0.], dist=0.)
        wz = Wall(normal=[0., 0., 1.], dist=0.)

        # (1)
        constraint_wy = espressomd.constraints.ShapeBasedConstraint(
            shape=wy, particle_type=1)
        S.constraints.add(constraint_wy)

        # (3)
        S.constraints.add(shape=wz, particle_type=2)

        # Check forces
        f_part = S.part[0].f

        self.assertEqual(f_part[0], 0.)
        self.assertEqual(f_part[1], 0.)
        self.assertEqual(f_part[2], 0.)

        S.integrator.run(0)
        f_part = S.part[0].f

        self.assertEqual(f_part[0], 0.)
        self.assertAlmostEqual(f_part[1], self.lj_force(
            eps=1.0, sig=1.0, r=1.21), places=10)
        self.assertAlmostEqual(f_part[2], self.lj_force(
            eps=1.5, sig=1.0, r=0.83), places=10)

        # Check removal
        for c in S.constraints:
            S.constraints.remove(c)

        for c in S.constraints:
            self.assertTrue(False)

if __name__ == "__main__":
    ut.main()
