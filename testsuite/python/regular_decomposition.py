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
import numpy as np
import itertools

np.random.seed(42)


class RegularDecomposition(ut.TestCase):
    system = espressomd.System(box_l=3 * [50.0])
    original_node_grid = tuple(system.cell_system.node_grid)

    def setUp(self):
        self.system.cell_system.set_regular_decomposition(
            use_verlet_lists=False)
        self.system.cell_system.node_grid = self.original_node_grid
        self.system.time_step = 1e-3

    def tearDown(self):
        self.system.part.clear()

    def check_resort(self):
        n_part = 2351

        # Add the particles on node 0, so that they have to be resorted
        particles = self.system.part.add(
            pos=n_part * [(0, 0, 0)], type=n_part * [1])

        # And now change their positions
        particles.pos = self.system.box_l * np.random.random((n_part, 3))

        # All particles should still be on node 0
        np.testing.assert_array_equal(np.copy(particles.node), 0)

        # Add an interacting particle in a corner of the box
        self.system.part.add(pos=(0.01, 0.01, 0.01), type=0)
        if espressomd.has_features(['LENNARD_JONES']):
            self.system.non_bonded_inter[0, 1].lennard_jones.set_params(
                epsilon=1.0, sigma=3.0, cutoff=6.0, shift=0.1)
            ref_energy = self.system.analysis.energy()['total']
            assert ref_energy > 10.

        # Distribute the particles on the nodes
        part_dist = self.system.cell_system.resort()

        # Check that we did not lose particles
        self.assertEqual(sum(part_dist), n_part + 1)

        # Check that we can still access all the particles
        # This basically checks if part_node and local_particles
        # are still in a valid state after the particle exchange
        self.assertEqual(sum(self.system.part.all().type), n_part)

        # Check that the system is still valid
        if espressomd.has_features(['LENNARD_JONES']):
            # energy calculation
            new_energy = self.system.analysis.energy()['total']
            self.assertEqual(new_energy, ref_energy)
        # force calculation
        self.system.integrator.run(0, recalc_forces=True)

        # Check particle transfer back to node 0
        old_nodes = np.copy(particles.node)
        particles.pos = n_part * [(0., 0., 0.)]
        new_nodes = np.copy(particles.node)
        np.testing.assert_array_equal(new_nodes, old_nodes)
        self.system.cell_system.resort()
        new_nodes = np.copy(particles.node)
        np.testing.assert_array_equal(new_nodes, 0)

    def test_resort(self):
        self.check_resort()

    @ut.skipIf(system.cell_system.get_state()["n_nodes"] != 4,
               "Skipping test: only runs for n_nodes == 4")
    def test_resort_alternating(self):
        # check particle resorting when the left and right cells are different
        self.system.cell_system.node_grid = [4, 1, 1]
        self.check_resort()

    def test_position_rounding(self):
        """This places a particle on the box boundary,
           with parameters that could cause problems with
           rounding."""
        self.system.box_l = [50.0, 50.0, 50.0]
        self.system.cell_system.skin = 0.4
        self.system.min_global_cut = 12.0 / 4.25
        self.system.part.add(pos=[25, 25, 0])
        self.assertEqual(1, len(self.system.part))

    @utx.skipIfMissingFeatures("LENNARD_JONES")
    def test_fully_connected_boundary(self):
        system = self.system
        system.part.clear()
        if system.cell_system.node_grid[1] != 1:
            ng = system.cell_system.node_grid
            system.cell_system.node_grid = [ng[0], 1, ng[2] * ng[1]]
        system.periodic = [True] * 3
        # Check that it's initially disabled
        self.assertEqual(system.cell_system.get_params()[
                         "fully_connected_boundary"], None)

        # check setting and getting the parameter
        system.cell_system.set_regular_decomposition(
            fully_connected_boundary=dict(direction="y", boundary="z"))
        self.assertEqual(system.cell_system.get_params()[
                         "fully_connected_boundary"], dict(direction="y", boundary="z"))
        # Check that the setting survives cell system re-initialization
        system.cell_system.min_global_cut = system.box_l / 4.1
        self.assertEqual(system.cell_system.get_params()[
                         "fully_connected_boundary"], dict(direction="y", boundary="z"))

        # Check particle visibility.
        # Place particles on a cubic lattice and use the
        # non_bonded_loop_trace() to check that all pairs are seen as expected
        fc_normal = np.array((0, 0, 1))  # z
        fc_normal_coord = 2  # z
        fc_dir = np.array((0, 1, 0))  # y
        N = 10
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            sigma=1, epsilon=1, cutoff=system.box_l[0] / N + 0.01, shift="auto")
        indices = [np.array((i, j, k)) for i in range(N)
                   for j in range(N) for k in range(N)]

        def id_for_idx(idx): return (
            idx[0] % N) * N * N + (idx[1] % N) * N + idx[2] % N

        ids = [id_for_idx(idx) for idx in indices]
        dx = system.box_l / N
        positions = [idx * dx for idx in indices]
        system.part.add(id=ids, pos=positions)
        particles = {i: system.part.by_id(i) for i in ids}

        def distance(id1, id2):
            return system.distance(
                particles[id1], particles[id2])
        distances = {tuple(i): distance(*i)
                     for i in itertools.combinations(ids, 2)}

        max_range = np.amax(system.box_l) / N
        two_cells = 2 * np.amax(system.cell_system.get_state()["cell_size"])
        two_cells_2d = two_cells * np.sqrt(2)
        two_cells_3d = two_cells * np.sqrt(3)
        assert np.all(system.box_l / 2 > two_cells)

        # next neighbors
        must_find_nn = [i for i, d in distances.items() if d <= max_range]

        # Fully connected neighbors
        indices_lower_boundary = [
            idx for idx in indices if idx[fc_normal_coord] == 0]
        must_find_fc = [tuple(sorted((id_for_idx(idx), id_for_idx(idx + i * fc_dir - fc_normal))))
                        for idx in indices_lower_boundary for i in range(-N + 1, N)]

        # all neighbors that must be found
        must_find = set(must_find_nn + must_find_fc)

        def assert_can_find(pair):
            # are the particles within a range that MAY be found by the
            # pair loop
            p1 = particles[pair[0]]
            p2 = particles[pair[1]]
            d = system.distance_vec(p1, p2)
            # if not accross periodic boundary: particles must be in cells
            # sharing at least one corner
            if np.abs(
                    p1.pos - p2.pos)[fc_normal_coord] < system.box_l[fc_normal_coord] / 2:
                self.assertLess(np.linalg.norm(d), two_cells_3d)
            # If across a the fully connected boundary
            # substract the distance in the fully connected direciont (all are
            # valid
            d_trans = d - d * fc_dir
            # in the other TWO directions, cells have to share a corner
            self.assertLess(np.linalg.norm(d_trans), two_cells_2d)

        # Use the cell system trace to get all pairs
        # as opposed to get_pairs() this does not have a distance check
        cs_pairs = system.cell_system.non_bonded_loop_trace()
        found = []
        for id1, id2, _rest1, _rest2, _rest3, _rest4 in cs_pairs:
            p = tuple(sorted((id1, id2)))  # Make the pair unique
            found.append(p)  # to check for double counting
            if p in must_find:
                must_find.remove(p)
            else:
                assert_can_find(p)  # close enough so that cells share a corner

        # Check for double counting of pairs
        self.assertEqual(len(found), len(set(found)))

        # check that all required pairs have been seen
        self.assertEqual(must_find, set([]))


if __name__ == "__main__":
    ut.main()
