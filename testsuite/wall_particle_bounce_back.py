from __future__ import division, print_function

import unittest as ut
import numpy

import espressomd
import espressomd.interactions
import espressomd.shapes as shapes
import tests_common


@ut.skipIf(not espressomd.has_features(["CONSTRAINTS", "LENNARD_JONES"]),
           "Features not available, skipping test!")
class ShapeBasedConstraintTest(ut.TestCase):

    box_l = 30.
    system = espressomd.System(box_l=3 * [box_l])

    def tearDown(self):
        self.system.part.clear()
        self.system.constraints.clear()

    def test_bounce_back(self):
        """Checks that a particle is bounced back from a wall it is heading towards.

        """
        system = self.system
        system.time_step = 0.01
        system.cell_system.skin = 0.4
        system.periodicity=[1,1,0] #make system non periodic in order to make bounce back boundary conditions work (it useses folded positions). Otherwise if the wall is placed on the periodic boundary the particle positions the bounce back method works with are folded back and then the existing code does not work. The code is also only guaranteed to work with flat walls (see implementation).
        system.part.add(pos=[0., 0., self.box_l-0.2], type=0, v=[0,0,1.0])
        
        system.non_bonded_inter[0, 1].lennard_jones.set_params(epsilon=0, sigma=0.0, cutoff=0.0, shift=0) #there needs to be an interaction for bounce back to work in the current implementation

        ceil = shapes.Wall(normal=[0, 0, -1], dist=-1.0*self.box_l)
        reflection_type=1
        c2 = system.constraints.add(particle_type=1, penetrable=int(False), only_positive=int(False), shape=ceil, reflection_type=reflection_type)

        const1 = system.constraints.add(c2)
        system.integrator.run(10)
        zpos=system.part[0].pos[2]
        self.assertAlmostEqual(
    zpos-29.9, 0.00, places=4,
                       msg="difference to analytical expected position too big")

if __name__ == "__main__":
    ut.main()
