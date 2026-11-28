#
# Copyright (C) 2010-2026 The ESPResSo project
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
import numpy as np
import espressomd


system = espressomd.System(box_l=[10., 10., 10.])
system.time_step = 0.01
system.cell_system.skin = 0.4
n_ranks = int(np.prod(system.cell_system.node_grid))
if n_ranks == 4:
    system.cell_system.node_grid = [4, 1, 1]


@ut.skipIf(n_ranks != 4,
           "Needs 4 MPI ranks: with fewer, every rank neighbors all others")
class ResortError(ut.TestCase):
    """
    A particle that leaves its MPI domain by more than one local box length
    must raise an error on all ranks. In particular, a failed particle sorting
    check must not leave the other ranks blocked in a collective call (the
    test would then run into the ctest timeout).
    """
    system = system

    def test_particle_leaves_domain(self):
        # local box length is 2.5: the particle ends up two domains away,
        # which is not reachable from a neighbor rank
        self.system.part.add(pos=[1., 1., 1.], v=[500., 0., 0.])
        with self.assertRaisesRegex(Exception, "moved more than one local box length"):
            self.system.integrator.run(1)


if __name__ == "__main__":
    ut.main()
