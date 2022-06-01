# Copyright (C) 2022 The ESPResSo project
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
import unittest as ut
import numpy as np
import espressomd
import espressomd.lees_edwards


class Test(ut.TestCase):
    system = espressomd.System(box_l=[1., 1., 1.])
    system.cell_system.skin = 0.01
    system.time_step = 0.01
    n_nodes = system.cell_system.get_state()["n_nodes"]
    node_grid = system.cell_system.node_grid

    def setUp(self):
        self.system.cell_system.set_regular_decomposition()
        self.system.cell_system.node_grid = self.node_grid
        self.system.periodicity = [True, True, True]

    def tearDown(self):
        self.system.part.clear()
        self.system.lees_edwards.protocol = None

    def get_neighbors_list(self, p1, dist):
        result = []
        for p2 in self.system.part.all():
            if p1.id != p2.id and self.system.distance(p1, p2) <= dist:
                result.append(p2.id)
        return result

    def check_neighbors(self, dist):
        for p in self.system.part.all():
            list_ref = np.sort(self.get_neighbors_list(p, dist))
            list_core = np.sort(self.system.cell_system.get_neighbors(p, dist))
            np.testing.assert_array_equal(list_core, list_ref)

    @ut.skipIf(n_nodes == 1, "only runs for 2 or more MPI ranks")
    def test_get_neighbors_different_domains(self):
        system = self.system
        get_neighbors = system.cell_system.get_neighbors
        pos_red = system.box_l - 0.01
        pos_black = np.array(3 * [0.01])
        with self.subTest(msg='no ghost update required'):
            # check particles on neighboring nodes
            p1, p2 = system.part.add(id=[1, 2], pos=[pos_black, pos_red])
            assert p1.node != p2.node
            np.testing.assert_array_equal(get_neighbors(p1, 0.1), [2])
            np.testing.assert_array_equal(get_neighbors(p2, 0.1), [1])
        system.part.clear()
        with self.subTest(msg='ghost update required'):
            # check particles on the same node
            p1, p2 = system.part.add(id=[1, 2], pos=[pos_black, pos_black])
            assert p1.node == p2.node
            np.testing.assert_array_equal(get_neighbors(p1, 0.1), [2])
            np.testing.assert_array_equal(get_neighbors(p2, 0.1), [1])
            # check after move to the same node
            p2.pos = pos_black + [0.12, 0., 0.]
            assert p1.node == p2.node
            np.testing.assert_array_equal(get_neighbors(p1, 0.1), [])
            np.testing.assert_array_equal(get_neighbors(p2, 0.1), [])
            # check after move to a neighboring node (ghost update required)
            p2.pos = pos_red
            assert p1.node == p2.node  # ghosts not updated yet
            np.testing.assert_array_equal(get_neighbors(p1, 0.1), [2])
            np.testing.assert_array_equal(get_neighbors(p2, 0.1), [1])
            # check after move to a neighboring node + forced ghost update
            system.integrator.run(0, recalc_forces=False)
            assert p1.node != p2.node
            np.testing.assert_array_equal(get_neighbors(p1, 0.1), [2])
            np.testing.assert_array_equal(get_neighbors(p2, 0.1), [1])

    @ut.skipIf(n_nodes == 3, "Only runs if number of MPI cores is not 3")
    def test_get_neighbors_random_positions(self):
        system = self.system
        system.part.add(pos=np.random.random((20, 3)) * system.box_l)
        for _ in range(10):
            system.part.all().pos = np.random.random((20, 3)) * system.box_l
            self.check_neighbors(0.2)
        system.part.clear()

    def test_n_square(self):
        """
        Check the N-square system is supported. The search distance is not
        bounded by the box size. No double counting occurs.
        """
        system = self.system
        system.cell_system.set_n_square()
        p1 = system.part.add(pos=[0., 0., 0.])
        p2 = system.part.add(pos=[0., 0., 0.49])
        self.assertEqual(system.cell_system.get_neighbors(p1, 0.5), [p2.id])
        self.assertEqual(system.cell_system.get_neighbors(p1, 50.), [p2.id])

    @ut.skipIf(n_nodes == 1, "Only runs if number of MPI cores is 2 or more")
    def test_domain_decomposition(self):
        """
        Check the neighbor search distance is bounded by the regular
        decomposition maximal range and that aperiodic systems are supported.
        """
        system = self.system
        get_neighbors = system.cell_system.get_neighbors
        local_box_l = np.copy(system.box_l / system.cell_system.node_grid)
        max_range = np.min(local_box_l) / 2.
        pos = np.array([0.01, 0.01, 0.01])
        p1, p2 = system.part.add(pos=[pos, -pos])
        self.check_neighbors(0.9999 * max_range)
        with np.testing.assert_raises_regex(ValueError, "pair search distance .* bigger than the decomposition range"):
            get_neighbors(p1, 1.0001 * max_range)

        # change periodicity
        pair_dist = np.linalg.norm(system.distance_vec(p1, p2))
        self.assertEqual(get_neighbors(p1, 1.1 * pair_dist), [p2.id])
        system.periodicity = [False, True, True]
        self.assertEqual(get_neighbors(p1, 1.1 * pair_dist), [])

    def test_hybrid_decomposition(self):
        """
        Set up colloids in a N-square cell system coupled to ions on a
        regular decomposition cell system.
        """
        system = self.system
        system.cell_system.set_hybrid_decomposition(
            n_square_types={1}, cutoff_regular=0.1)
        p_ion = system.part.add(pos=[0., 0., 0.], type=0)
        p_colloid = system.part.add(pos=[0., 0., 0.], type=1)

        msg = "Cannot search for neighbors in the hybrid decomposition cell system"
        with np.testing.assert_raises_regex(RuntimeError, msg):
            system.cell_system.get_neighbors(p_ion, 0.05)
        with np.testing.assert_raises_regex(RuntimeError, msg):
            system.cell_system.get_neighbors(p_colloid, 0.05)

    def test_lees_edwards(self):
        """
        Check the Lees-Edwards position offset is taken into account
        in the distance calculation.
        """
        system = self.system
        protocol = espressomd.lees_edwards.LinearShear(
            initial_pos_offset=0.1, time_0=0., shear_velocity=1.0)
        system.lees_edwards.set_boundary_conditions(
            shear_direction=0, shear_plane_normal=2, protocol=protocol)
        p1 = system.part.add(pos=[0., 0., 0.])
        p2 = system.part.add(pos=[0., 0., 0.99])
        self.assertEqual(system.cell_system.get_neighbors(p1, 0.10), [])
        self.assertEqual(system.cell_system.get_neighbors(p1, 0.15), [p2.id])


if __name__ == "__main__":
    ut.main()
