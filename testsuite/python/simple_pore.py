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
import unittest as ut
import unittest_decorators as utx
import espressomd
from espressomd.shapes import SimplePore, Cylinder

# Integration test for simple pore
# The rationale is to hit the pore everywhere with particles
# and check that it does not blow up. The cylinder is needed
# because the pore is tilted with respect to the box, without
# it particles could enter the constraint over the periodic boundaries,
# leading to force jumps.


@utx.skipIfMissingFeatures(["LENNARD_JONES"])
class SimplePoreConstraint(ut.TestCase):

    def test_orientation(self):
        pore = SimplePore(axis=[1., 0., 0.], radius=2., smoothing_radius=.1,
                          length=2., center=[5., 5., 5.])

        d, _ = pore.calc_distance(position=[.0, .0, .0])
        self.assertGreater(d, 0.)

        d, _ = pore.calc_distance(position=[5., 5., .0])
        self.assertLess(d, 0.)

    def test_stability(self):
        s = espressomd.System(box_l=[1.0, 1.0, 1.0])
        box_yz = 15.
        box_x = 20.
        s.box_l = [box_x, box_yz, box_yz]
        s.time_step = 0.01
        s.cell_system.skin = 0.4

        lj_eps = 1.0
        lj_sig = 1.0
        lj_cut = lj_sig * 2**(1. / 6.)

        s.constraints.add(
            particle_type=0, penetrable=False, only_positive=False,
            shape=SimplePore(
                axis=[1., 0.5, 0.5], radius=3., smoothing_radius=.1,
                length=5, center=[.5 * box_x, .5 * box_yz, .5 * box_yz]))
        s.constraints.add(
            particle_type=0, penetrable=False, only_positive=False,
            shape=Cylinder(
                axis=[1., 0, 0], radius=0.5 * box_yz, length=4 * lj_cut + box_x,
                center=[.5 * box_x, .5 * box_yz, .5 * box_yz], direction=-1))

        s.non_bonded_inter[0, 1].lennard_jones.set_params(
            epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")

        for i in range(200):
            rpos = [i * (box_x / 200.), 0.5 * box_yz, 0.5 * box_yz]
            s.part.add(id=i, pos=rpos, type=1, v=[1., 1., 1.])

        start_energy = s.analysis.energy()['total']
        s.integrator.run(1000)
        end_energy = s.analysis.energy()['total']
        rel_diff = abs(end_energy - start_energy) / start_energy

        self.assertLess(rel_diff, 1e-3)


if __name__ == "__main__":
    ut.main()
