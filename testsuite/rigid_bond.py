
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
# Tests particle property setters/getters
from __future__ import print_function
import unittest as ut
import espressomd
import numpy as np
from espressomd.interactions import RigidBond


@ut.skipIf(not espressomd.has_features("BOND_CONSTRAINT"),
           "Test requires BOND_CONSTRAINT feature")
class RigidBondTest(ut.TestCase):

    def test(self):
        target_acc = 1E-3
        tol = 1.2 * target_acc
        s = espressomd.System(box_l=[1.0, 1.0, 1.0])
        s.seed = s.cell_system.get_state()['n_nodes'] * [1234]
        s.box_l = 10, 10, 10
        s.cell_system.skin = 0.4
        s.time_step = 0.01
        s.thermostat.set_langevin(kT=1, gamma=1)
        r = RigidBond(r=1.2, ptol=1E-3, vtol=target_acc)
        s.bonded_inter.add(r)

        for i in range(5):
            s.part.add(id=i, pos=(i * 1.2, 0, 0))
            if i > 0:
                s.part[i].bonds = ((r, i - 1),)
        s.integrator.run(5000)
        for i in range(1, 5):
            d = s.distance(s.part[i], s.part[i - 1])
            v_d = s.distance_vec(s.part[i], s.part[i - 1])
            self.assertLess(abs(d - 1.2), tol)
            # Velocity projection on distance vector
            vel_proj = np.dot(s.part[i].v - s.part[i - 1].v, v_d) / d
            self.assertLess(vel_proj, tol)


if __name__ == "__main__":
    #print("Features: ", espressomd.features())
    ut.main()
