#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
# Tests hat potential
from __future__ import print_function
import unittest as ut
import espressomd
import numpy as np

@ut.skipIf(not espressomd.has_features("HAT"),"Skipped because feature is disabled")
class HatTest(ut.TestCase):
    def force(self, F_max, r_cut, r):
        if r > 0 and r < r_cut:
            return F_max * (1. - r/r_cut )
        else:
            return 0.

    def pot(self, F_max, r_cut, r):
        if r < r_cut:
            return F_max * (r - r_cut) * ( (r + r_cut) / (2. * r_cut) - 1.)
        else:
            return 0.

    def test(self):
        s = espressomd.System(box_l=[1.0, 1.0, 1.0])
        s.seed = s.cell_system.get_state()['n_nodes'] * [1234]
        s.box_l = 3 * [10]
        s.time_step = 0.01
        s.cell_system.skin=0.4

        s.part.add(id=0, pos=0.5 * s.box_l, type=0)
        s.part.add(id=1, pos=0.5 * s.box_l, type=0)

        F_max =3.145
        cutoff = 1.3

        s.non_bonded_inter[0,0].hat.set_params(F_max=F_max, cutoff=cutoff)

        dx = cutoff / 90.
        r0 = 0.5 * s.box_l[0]

        for i in range(100):
            r = r0 - i*dx
            s.part[1].pos=[r, 0.5 * s.box_l[1], 0.5*s.box_l[2]]
            s.integrator.run(0)
            self.assertAlmostEqual(self.force(F_max, cutoff, i*dx), s.part[0].f[0], places=7)
            self.assertAlmostEqual(self.pot(F_max, cutoff, i*dx), s.analysis.energy()['total'], places=7)

if __name__ == "__main__":
    ut.main()
