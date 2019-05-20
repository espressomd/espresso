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
if espressomd.has_features("LB_WALBERLA"):
    from espressomd.lb import LBFluidWalberla
import numpy as np


@ut.skipIf(not espressomd.has_features("LB_WALBERLA"), "Skipping for LACK of LB_WALBERLA")
class LbWalberlaTest(ut.TestCase):

    def test(self):
        s = espressomd.System(box_l=(9.6, 12, 15.6))
        s.time_step = 0.2
        s.cell_system.skin = 0.
        dens_init = 1.3
        lbf = LBFluidWalberla(
            agrid=.6,
            dens=dens_init,
            visc=2.5,
            tau=s.time_step)
        s.actors.add(lbf)
        max_ind = s.box_l / .6
        for i in range(int(max_ind[0])):
            for j in range(int(max_ind[1])):
                for k in range(int(max_ind[2])):
                    assert np.linalg.norm(lbf[i, j, k].velocity) < 1E-15
                    v = np.array((i * .5, j * 1.5, k * 2.5))
                    lbf[i, j, k].velocity = v
                    assert np.linalg.norm(lbf[i, j, k].velocity - v) < 1E-10

                    assert abs(
                        lbf[i, j, k].density - dens_init / 0.6**3) < 1E-10
                    rho = i * j * k * 0.5 + 0.8
                    lbf[i, j, k].density = rho
                    assert abs(lbf[i, j, k].density - rho) < 1E-10

                    pop = np.array((i * j * k, i, -i, j, -j, k, -k,
                                    i + j, i - j, -i + j, -i - j, i + k, i - k, -i + k, -i - k, j + k, j - k, -j + k, -j - k))
                    lbf[i, j, k].population = pop
                    lb_pop = lbf[i, j, k].population
                    for p_idx in range(len(pop)):
                        assert abs(lb_pop[p_idx] - pop[p_idx]) < 1E-10

        s.actors.remove(lbf)


if __name__ == "__main__":
    # print("Features: ", espressomd.features())
    ut.main()
