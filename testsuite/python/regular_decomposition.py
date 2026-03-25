#
# Copyright (C) 2013-2026 The ESPResSo project
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

        def id_for_idx(idx):
            return int((idx[0] % N) * N * N + (idx[1] % N) * N + idx[2] % N)

        ids = [id_for_idx(idx) for idx in indices]
        dx = system.box_l / N
        positions = [idx * dx for idx in indices]
        system.part.add(id=ids, pos=positions)
        particles = {i: system.part.by_id(i) for i in ids}
        pos_folded = {pid: p.pos_folded for pid, p in particles.items()}

        vec_matrix = {
            (pid1, pid2): system.distance_vec(pos_folded[pid1], pos_folded[pid2])
            for pid1, pid2 in itertools.combinations(ids, 2)}
        dist_sq_matrix = {k: np.dot(v, v) for k, v in vec_matrix.items()}

        max_range = np.amax(system.box_l) / N
        max_range_sq = max_range**2
        two_cells = 2 * np.amax(system.cell_system.get_state()["cell_size"])
        two_cells_2d = two_cells * np.sqrt(2)
        two_cells_3d = two_cells * np.sqrt(3)
        assert np.all(system.box_l / 2 > two_cells)

        # next neighbors
        must_find_nn = [
            pid for pid, dist_sq in dist_sq_matrix.items() if dist_sq <= max_range_sq]

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
            # if not across periodic boundary: particles must be in cells
            # sharing at least one corner
            d_abs = np.abs(p1.pos - p2.pos)
            if d_abs[fc_normal_coord] < system.box_l[fc_normal_coord] / 2:
                self.assertLess(dist_sq_matrix[pair], two_cells_3d**2)
            # If across a fully connected boundary, subtract the distance
            # in the fully connected direction (all are valid)
            d = vec_matrix[pair]
            d_trans = d - d * fc_dir
            # in the other TWO directions, cells have to share a corner
            self.assertLess(np.dot(d_trans, d_trans), two_cells_2d**2)

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

    @utx.skipIfMissingFeatures("LENNARD_JONES")
    def test_fully_connected_boundary_two_particles(self):
        """Place two particles on opposite sides of the fully connected
        normal direction (z) with random positions along the fully connected
        direction (y) and x=0. Check that the pair is always discovered
        by the non-bonded loop trace."""
        system = self.system
        system.part.clear()
        if system.cell_system.node_grid[1] != 1:
            ng = system.cell_system.node_grid
            system.cell_system.node_grid = [ng[0], 1, ng[2] * ng[1]]

        # Use an asymmetric box: box_l[y] > box_l[z] so that there are
        # more cells along the fully connected direction (y) than along
        # the boundary normal (z).
        system.box_l = [50.0, 50.0, 30.0]
        system.cell_system.set_regular_decomposition(
            fully_connected_boundary=dict(direction="y", boundary="z"))

        skin = 0.4
        system.cell_system.skin = skin
        cutoff = 2 * skin + 0.01
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            sigma=1, epsilon=1, cutoff=cutoff, shift="auto")

        rng = np.random.RandomState(seed=1234)

        for i in range(1000):
            system.part.clear()
            # Random y positions (fully connected direction)
            y1 = rng.uniform(0, system.box_l[1])
            y2 = rng.uniform(0, system.box_l[1])
            # Opposite sides of z (fully connected normal direction)
            z1 = rng.uniform(0, skin)
            z2 = rng.uniform(system.box_l[2] - skin, system.box_l[2])
            # x = 0 (3rd direction)
            system.part.add(id=0, pos=[0, y1, z1])
            system.part.add(id=1, pos=[0, y2, z2])

            cs_pairs = system.cell_system.non_bonded_loop_trace()
            found = set()
            for id1, id2, *_rest in cs_pairs:
                found.add(tuple(sorted((id1, id2))))
            self.assertIn(
                (0, 1), found,
                msg=f"Pair not found for positions {[0, y1, z1]} "
                f"and {[0, y2, z2]} (iteration {i})")

    @utx.skipIfMissingFeatures("LENNARD_JONES")
    def test_fully_connected_boundary_at_boundary(self):
        """Place two particles on opposite sides of the fully connected
        normal direction (y) with random positions along the fully connected
        direction (z) and x=0. Uses direction='z', boundary='y' so that
        the linear index ordering exposes the at_boundary check bug:
        at_boundary uses global_size[coord] (a ghost index) instead of
        global_size[coord]-1 (the last real cell), so only the lower
        boundary gets fully connected neighbors.

        The column-major linear index x + sx*(y + sy*z) gives z the
        highest weight. When the y=0 particle has lower z, it has the
        lower linear index and discovers the pair via its red (higher
        index) neighbors where FC is active. When the y=box_l particle
        has the lower z, the pair must be discovered from the y=box_l
        side, which requires the upper boundary to also have FC."""
        system = self.system
        system.part.clear()
        system.box_l = [50.0, 50.0, 50.0]
        if system.cell_system.node_grid[2] != 1:
            ng = system.cell_system.node_grid
            system.cell_system.node_grid = [ng[0] * ng[2], ng[1], 1]

        system.cell_system.set_regular_decomposition(
            fully_connected_boundary=dict(direction="z", boundary="y"))

        skin = 0.4
        system.cell_system.skin = skin
        cutoff = 2 * skin + 0.01
        system.non_bonded_inter[0, 0].lennard_jones.set_params(
            sigma=1, epsilon=1, cutoff=cutoff, shift="auto")

        def get_pairs():
            cs_pairs = system.cell_system.non_bonded_loop_trace()
            found = set()
            for id1, id2, *_rest in cs_pairs:
                found.add(tuple(sorted((id1, id2))))
            return found

        rng = np.random.RandomState(seed=1234)

        for i in range(1000):
            z1 = rng.uniform(0, system.box_l[2])
            z2 = rng.uniform(0, system.box_l[2])
            y_bot = rng.uniform(0, skin)
            y_top = rng.uniform(system.box_l[1] - skin, system.box_l[1])
            z_lo, z_hi = min(z1, z2), max(z1, z2)

            # Check bottom boundary FC: y=0 particle at lower z gives it
            # the lower linear index, so the pair is discovered from the
            # y=0 (bottom) cell's red neighbor list.
            system.part.clear()
            system.part.add(id=0, pos=[0, y_bot, z_lo])
            system.part.add(id=1, pos=[0, y_top, z_hi])
            self.assertIn(
                (0, 1), get_pairs(),
                msg=f"Bottom boundary: pair not found for "
                f"{[0, y_bot, z_lo]} and {[0, y_top, z_hi]} "
                f"(iteration {i})")

            # Check top boundary FC: y=box_l particle at lower z gives it
            # the lower linear index, so the pair is discovered from the
            # y=box_l (top) cell's red neighbor list.
            system.part.clear()
            system.part.add(id=0, pos=[0, y_bot, z_hi])
            system.part.add(id=1, pos=[0, y_top, z_lo])
            self.assertIn(
                (0, 1), get_pairs(),
                msg=f"Top boundary: pair not found for "
                f"{[0, y_bot, z_hi]} and {[0, y_top, z_lo]} "
                f"(iteration {i})")

    def test_fully_connected_boundary_periodicity_check(self):
        """Check that setting up a fully connected boundary raises an error
        when the boundary normal direction is non-periodic."""
        system = self.system
        system.part.clear()
        old_periodicity = list(system.periodicity)
        # Ensure node_grid is compatible with fc direction="y"
        # (y must have 1 rank) so the node_grid check passes first
        if system.cell_system.node_grid[1] != 1:
            ng = system.cell_system.node_grid
            system.cell_system.node_grid = [ng[0], 1, ng[1] * ng[2]]
        try:
            system.periodicity = [True, True, False]
            with self.assertRaisesRegex(
                    RuntimeError,
                    "fully connected boundary requires periodicity"):
                system.cell_system.set_regular_decomposition(
                    fully_connected_boundary=dict(
                        direction="y", boundary="z"))
        finally:
            system.periodicity = old_periodicity
            system.cell_system.set_regular_decomposition()


if __name__ == "__main__":
    ut.main()
