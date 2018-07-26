from __future__ import division, print_function
import unittest as ut
import espressomd
from espressomd.shapes import SimplePore,Cylinder

# Integration test for simple pore
# The radional is to hit the pore everywhere with particles
# and check that it does not blow up. The cylinder is needed
# because the pore is tilted with respect to the box, without
# it particles could enter the constraint over the periodic boundarys,
# leading to force jumps.

@ut.skipIf(not espressomd.has_features(["LENNARD_JONES"]),
           "Features not available, skipping test!")
class SimplePoreConstraint(ut.TestCase):
    def test(self):
        s = espressomd.System(box_l=[1.0, 1.0, 1.0])
        s.seed = s.cell_system.get_state()['n_nodes'] * [1234]
        box_yz = 15.
        box_x = 20.
        s.box_l = [box_x, box_yz, box_yz]
        s.time_step = 0.01
        s.cell_system.skin = 0.4

        lj_eps = 1.0
        lj_sig = 1.0
        lj_cut = lj_sig * 2**(1./6.)

        s.constraints.add(particle_type=0, penetrable=False, only_positive=False, shape=SimplePore(axis=[1.,0.5,0.5], radius=3., smoothing_radius=.1,length=5, center=[.5*box_x, .5*box_yz,.5*box_yz]))
        s.constraints.add(particle_type=0, penetrable=False, only_positive=False, shape=Cylinder(axis=[1.,0.0,0], radius=0.5*box_yz,length=4*lj_cut+box_x, center=[.5*box_x, .5*box_yz,.5*box_yz],direction=-1))

        s.non_bonded_inter[0, 1].lennard_jones.set_params(
            epsilon=lj_eps, sigma=lj_sig,
            cutoff=lj_cut, shift="auto")

        for i in range(200):
            rpos=[i*(box_x / 200.), 0.5 * box_yz, 0.5 * box_yz]
            s.part.add(id=i, pos=rpos, type=1, v=[1.,1.,1.])

        start_energy = s.analysis.energy()['total']
        s.integrator.run(1000)
        end_energy = s.analysis.energy()['total']
        rel_diff = abs(end_energy-start_energy) / start_energy

        self.assertLess(rel_diff, 1e-3)

if __name__ == "__main__":
    ut.main()
