#
# Copyright (C) 2013-2022 The ESPResSo project
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
import unittest as ut
import unittest_decorators as utx
import espressomd
import espressomd.interactions
import numpy as np
import itertools


@utx.skipIfMissingFeatures("BOND_CONSTRAINT")
class RigidBondTest(ut.TestCase):

    def test(self):
        target_acc = 1E-3
        tol = 1.2 * target_acc
        system = espressomd.System(box_l=[10., 10., 10.])
        system.cell_system.skin = 0.4
        system.time_step = 0.01
        system.thermostat.set_langevin(kT=1, gamma=1, seed=42)
        rigid_bond = espressomd.interactions.RigidBond(
            r=1.2, ptol=1E-3, vtol=target_acc)
        system.bonded_inter.add(rigid_bond)

        # create polymer
        last_p = None
        for i in range(5):
            p = system.part.add(pos=(i * 1.2, 0, 0))
            if last_p is not None:
                p.add_bond((rigid_bond, last_p))
            last_p = p

        system.integrator.run(5000)

        # check every bond
        p1_iter, p2_iter = itertools.tee(system.part)
        next(p2_iter, None)  # advance second iterator by 1 step
        for p1, p2 in zip(p1_iter, p2_iter):
            d = system.distance(p2, p1)
            v_d = system.distance_vec(p2, p1)
            self.assertAlmostEqual(d, 1.2, delta=tol)
            # Velocity projection on distance vector
            vel_proj = np.dot(p2.v - p1.v, v_d) / d
            self.assertLess(vel_proj, tol)


if __name__ == "__main__":
    ut.main()
