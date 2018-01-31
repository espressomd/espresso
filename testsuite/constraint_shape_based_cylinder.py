# Tests if shape based constraints can be added to a system both by
#  (1) defining a constraint object which is then added
#  (2) and via keyword arguments.
# Checks, that cylinder constraints with LJ interactions exert forces
# on a test particle (that is, the constraints do what they should).

from __future__ import division, print_function

import unittest as ut

import espressomd

from espressomd import interactions
import espressomd.shapes as shapes
import tests_common


@ut.skipIf(not espressomd.has_features(["CONSTRAINTS", "LENNARD_JONES"]),
           "Features not available, skipping test!")
class ShapeBasedConstraintTest(ut.TestCase):
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

    def test_cylinder(self):
        system = espressomd.System()
        self.prepare(system)

        # check force calculation of cylinder constraint
        system.part[0].pos = [self.box_l / 2.0, 1.02, self.box_l / 2.0]
        interactio_dir = -1  # constraint is directed inwards
        cylinder_shape = shapes.Cylinder(center=[self.box_l / 2.0, self.box_l / 2.0, self.box_l / 2.0], axis=[
                                         0, 0, 1], direction=interactio_dir, radius=self.box_l / 2.0, length=self.box_l + 5)  # +5 in order to have no top or bottom
        penetrability = 0  # inpenetrable
        outer_cylinder_constraint = espressomd.constraints.ShapeBasedConstraint(
            shape=cylinder_shape, particle_type=1, penetrable=penetrability)
        outer_cylinder_wall = system.constraints.add(outer_cylinder_constraint)
        system.integrator.run(0)  # update forces

        # test forces on walls
        self.assertAlmostEqual(-1.0 * outer_cylinder_wall.total_force()[1], tests_common.lj_force(espressomd, cutoff=2.0, offset=0.,
                                                                                                  eps=1.0, sig=1.0, r=1.02), places=10)  # minus for newtons thrid law


if __name__ == "__main__":
    ut.main()
