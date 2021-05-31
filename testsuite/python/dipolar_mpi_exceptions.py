#
# Copyright (C) 2021 The ESPResSo project
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

import espressomd.magnetostatics


@utx.skipIfMissingFeatures(["DIPOLES"])
class MagnetostaticsInterface(ut.TestCase):
    system = espressomd.System(box_l=[10., 10., 10.])
    n_nodes = system.cell_system.get_state()["n_nodes"]

    def tearDown(self):
        self.system.actors.clear()

    @ut.skipIf(n_nodes == 1, "only runs for 2+ MPI ranks")
    def test_exceptions_mpi(self):
        actor = espressomd.magnetostatics.DipolarDirectSumCpu(prefactor=1)
        with self.assertRaisesRegex(RuntimeError, 'MPI parallelization not supported'):
            self.system.actors.add(actor)
        self.system.actors.clear()
        actor = espressomd.magnetostatics.DipolarDirectSumWithReplicaCpu(
            prefactor=1, n_replica=2)
        with self.assertRaisesRegex(RuntimeError, 'MPI parallelization not supported'):
            self.system.actors.add(actor)


if __name__ == "__main__":
    ut.main()
