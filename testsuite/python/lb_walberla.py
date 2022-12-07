# Copyright (C) 2013-2018 The ESPResSo project
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
if espressomd.has_features("WALBERLA"):
    from espressomd.lb import LBFluidWalberla
import numpy as np


@ut.skipIf(not espressomd.has_features("WALBERLA"),
           "Skipping for LACK of WALBERLA")
class LbWalberlaTest(ut.TestCase):

    def test(self):
        s = espressomd.System(box_l=(12 * .6, 18 * .6, 12 * .6))
        s.time_step = 0.2
        s.cell_system.skin = 0.
        dens_init = 1.3
        lbf = LBFluidWalberla(
            agrid=.6,
            density=dens_init,
            viscosity=2.5,
            tau=s.time_step)
        s.actors.add(lbf)
        max_ind = s.box_l / .6
        for i in range(int(max_ind[0])):
            for j in range(int(max_ind[1])):
                for k in range(int(max_ind[2])):
                    # velocity
                    np.testing.assert_allclose(
                        lbf[i, j, k].velocity, [0, 0, 0], atol=1E-15)
                    v = np.array((i * .5, j * 1.5, k * 2.5))
                    lbf[i, j, k].velocity = v
                    np.testing.assert_allclose(
                        lbf[i, j, k].velocity, v, atol=1E-10)
                    # density
                    self.assertAlmostEqual( 
                        lbf[i, j, k].density, dens_init, delta=1E-10)
                    rho = i * j * k * 0.5 + 0.8
                    lbf[i, j, k].density = rho
                    self.assertAlmostEqual(
                        lbf[i, j, k].density, rho, delta=1E-10)
                    # population
                    pop = np.array((i * j * k, i, -i, j, -j, k, -k,
                                    i + j, i - j, -i + j, -i - j, i + k,
                                    i - k, -i + k, -i - k, j + k, j - k,
                                    -j + k, -j - k), dtype=float)
                    lbf[i, j, k].population = pop
                    lb_pop = lbf[i, j, k].population
                    np.testing.assert_allclose(lb_pop, pop, atol=1E-10)
                    # last applied force
                    last_applied_force = np.array((i, -j, 1 + k))
                    lbf[i, j, k].last_applied_force = last_applied_force
                    lb_last_applied_force = lbf[i, j, k].last_applied_force
                    np.testing.assert_allclose(
                        lb_last_applied_force, last_applied_force, atol=1E-10)

        s.actors.remove(lbf)


if __name__ == "__main__":
    ut.main()
