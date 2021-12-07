#
# Copyright (C) 2013-2019 The ESPResSo project
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
import espressomd
import numpy as np

np.random.seed(42)


class NSquare(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])

    def setUp(self):
        self.system.cell_system.set_n_square(use_verlet_lists=False)
        self.system.time_step = 1e-3
        self.system.cell_system.skin = 0.15

    def tearDown(self):
        self.system.part.clear()

    def test_load_balancing(self):
        n_part = 235
        n_nodes = self.system.cell_system.get_state()['n_nodes']
        n_part_avg = n_part // n_nodes

        # Add the particles on node 0, so that they have to be resorted
        partcls = self.system.part.add(
            pos=n_part * [(0, 0, 0)], type=n_part * [1])

        # And now change their positions
        partcls.pos = self.system.box_l * \
            np.random.random((n_part, 3))

        # Add an interacting particle in a corner of the box
        self.system.part.add(pos=[(0.01, 0.01, 0.01)], type=[0])
        if espressomd.has_features(['LENNARD_JONES']):
            self.system.non_bonded_inter[0, 1].lennard_jones.set_params(
                epsilon=1.0, sigma=0.14, cutoff=0.15, shift=0.1)
            ref_energy = self.system.analysis.energy()['total']
            assert ref_energy > 10.

        # Distribute the particles on the nodes
        part_dist = self.system.cell_system.resort()

        # Check that we did not lose particles
        self.assertEqual(sum(part_dist), n_part + 1)

        # Check that the particles are evenly distributed
        for node_parts in part_dist:
            self.assertAlmostEqual(node_parts, n_part_avg, delta=2)

        # Check that we can still access all the particles
        # This basically checks if part_node and local_particles
        # are still in a valid state after the particle exchange
        self.assertEqual(sum(partcls.type), n_part)

        # Check that the system is still valid
        if espressomd.has_features(['LENNARD_JONES']):
            # energy calculation
            new_energy = self.system.analysis.energy()['total']
            self.assertEqual(new_energy, ref_energy)
        # force calculation
        self.system.integrator.run(0, recalc_forces=True)


if __name__ == "__main__":
    ut.main()
