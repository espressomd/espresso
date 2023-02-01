#
# Copyright (C) 2017-2022 The ESPResSo project
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

import numpy as np
import unittest as ut
import unittest_decorators as utx
import espressomd


@utx.skipIfMissingFeatures("EXCLUSIONS")
class Test(ut.TestCase):
    """
    Check the auto-exclusions feature against various polymer topologies.
    Verify the exclusion list is still correct when the polymer spans
    three contiguous, yet different, cells, such that particles at one
    end of the polymer don't have access to the particles at the other
    end of the polymer via the ghost layer.
    """
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    system.cell_system.skin = 0.1
    node_grid_ref = list(system.cell_system.node_grid)

    def tearDown(self):
        self.system.part.clear()
        self.system.cell_system.node_grid = self.node_grid_ref

    def set_particles_on_cube(self):
        # place particles on a cube centered at the origin
        length = self.system.cell_system.skin / 2.
        self.system.part.add(pos=length * (np.array([0, 0, 0]) - 1))
        self.system.part.add(pos=length * (np.array([2, 0, 0]) - 1))
        self.system.part.add(pos=length * (np.array([2, 2, 0]) - 1))
        self.system.part.add(pos=length * (np.array([0, 2, 0]) - 1))
        self.system.part.add(pos=length * (np.array([0, 2, 2]) - 1))
        self.system.part.add(pos=length * (np.array([2, 2, 2]) - 1))
        self.system.part.add(pos=length * (np.array([2, 0, 2]) - 1))
        self.system.part.add(pos=length * (np.array([0, 0, 2]) - 1))

        # particles should be equally distributed on all nodes
        n_nodes = np.prod(self.node_grid_ref)
        if n_nodes in (1, 2, 4, 8):
            p_nodes = sorted(list(self.system.part.all().node))
            p_nodes_ref = list(np.repeat(np.arange(n_nodes), 8 // n_nodes))
            assert p_nodes == p_nodes_ref

    def set_particles_on_line(self, length):
        box_l_z = self.system.box_l[2]
        for i in range(length):
            self.system.part.add(pos=np.array([0, 0, box_l_z]) * i / length)

        # particles should be distributed on multiple nodes
        n_nodes = np.prod(self.node_grid_ref)
        if n_nodes > 1:
            assert len(set(self.system.part.all().node)) > 1

    def test_linear(self):
        bond = espressomd.interactions.Virtual()
        system = self.system
        system.bonded_inter.add(bond)

        def check():
            # topology: 0-1-2-...-n
            length = len(system.part)
            for pid in range(length - 1):
                system.part.by_id(pid).add_bond((bond, pid + 1))
            system.auto_exclusions(distance=1)

            for pid in range(1, length - 1):
                excl = sorted(list(system.part.by_id(pid).exclusions))
                self.assertEqual(excl, [pid - 1, pid + 1])

            excl = list(system.part.by_id(0).exclusions)
            self.assertEqual(excl, [1])

            excl = list(system.part.by_id(length - 1).exclusions)
            self.assertEqual(excl, [length - 2])

        with self.subTest(msg='cube'):
            self.set_particles_on_cube()
            check()

        self.tearDown()

        with self.subTest(msg='line'):
            system.cell_system.node_grid = [1, 1, np.prod(self.node_grid_ref)]
            self.set_particles_on_line(16)
            check()

    def test_ring(self):
        bond = espressomd.interactions.Virtual()
        system = self.system
        system.bonded_inter.add(bond)

        def check():
            # topology: 0-1-2-...-n-0
            length = len(system.part)
            for pid in range(length):
                system.part.by_id(pid).add_bond((bond, (pid + 1) % length))
            system.auto_exclusions(distance=2)

            for pid in range(length):
                excl = sorted(list(system.part.by_id(pid).exclusions))
                excl_ref = np.mod([pid - 1, pid - 2, pid + 1, pid + 2], length)
                self.assertEqual(excl, sorted(list(excl_ref)))

        with self.subTest(msg='cube'):
            self.set_particles_on_cube()
            check()

        self.tearDown()

        with self.subTest(msg='line'):
            system.cell_system.node_grid = [1, 1, np.prod(self.node_grid_ref)]
            self.set_particles_on_line(16)
            check()

    def test_branched(self):
        bond = espressomd.interactions.Virtual()
        system = self.system
        system.bonded_inter.add(bond)

        length = system.cell_system.skin / 2.
        p0 = system.part.add(pos=length * (np.array([0, 0, 0]) - 1))
        p1 = system.part.add(pos=length * (np.array([2, 0, 0]) - 1))
        p2 = system.part.add(pos=length * (np.array([0, 2, 0]) - 1))
        p3 = system.part.add(pos=length * (np.array([2, 2, 0]) - 1))

        # topology: 0-1(-2)-3
        p1.add_bond((bond, p0))
        p2.add_bond((bond, p1))
        p3.add_bond((bond, p1))
        system.auto_exclusions(distance=2)

        self.assertEqual(sorted(list(p0.exclusions)), [1, 2, 3])
        self.assertEqual(sorted(list(p1.exclusions)), [0, 2, 3])
        self.assertEqual(sorted(list(p2.exclusions)), [0, 1, 3])
        self.assertEqual(sorted(list(p3.exclusions)), [0, 1, 2])

    def test_diamond(self):
        bond = espressomd.interactions.Virtual()
        system = self.system
        system.bonded_inter.add(bond)

        length = system.cell_system.skin / 2.
        p0 = system.part.add(pos=length * (np.array([0, 0, 2]) - 1))
        p1 = system.part.add(pos=length * (np.array([0, 0, 0]) - 1))
        p2 = system.part.add(pos=length * (np.array([2, 0, 0]) - 1))
        p3 = system.part.add(pos=length * (np.array([0, 2, 0]) - 1))
        p4 = system.part.add(pos=length * (np.array([2, 2, 0]) - 1))
        p5 = system.part.add(pos=length * (np.array([2, 2, 2]) - 1))

        # topology: 0-1(-2-4-5)-3-4-5
        p1.add_bond((bond, p0))
        p2.add_bond((bond, p1))
        p3.add_bond((bond, p1))
        p2.add_bond((bond, p4))
        p3.add_bond((bond, p4))
        p5.add_bond((bond, p4))
        system.auto_exclusions(distance=3)

        self.assertEqual(sorted(list(p0.exclusions)), [1, 2, 3, 4])
        self.assertEqual(sorted(list(p1.exclusions)), [0, 2, 3, 4, 5])
        self.assertEqual(sorted(list(p2.exclusions)), [0, 1, 3, 4, 5])
        self.assertEqual(sorted(list(p3.exclusions)), [0, 1, 2, 4, 5])
        self.assertEqual(sorted(list(p4.exclusions)), [0, 1, 2, 3, 5])
        self.assertEqual(sorted(list(p5.exclusions)), [1, 2, 3, 4])

    def test_id_gaps(self):
        bond = espressomd.interactions.Virtual()
        system = self.system
        system.bonded_inter.add(bond)

        p0 = system.part.add(id=0, pos=[0.1, 0.1, 0.1])
        p4 = system.part.add(id=4, pos=[0.2, 0.1, 0.1])

        # topology: 0-4
        p0.add_bond((bond, p4))
        system.auto_exclusions(distance=1)

        self.assertEqual(list(p0.exclusions), [4])
        self.assertEqual(list(p4.exclusions), [0])

    def test_non_pairwise(self):
        """
        Check that bonds with 0, 2 and 3 partners don't generate exclusions.
        """
        angle = espressomd.interactions.AngleHarmonic(bend=1., phi0=2.)
        dihe = espressomd.interactions.Dihedral(bend=1., mult=2, phase=2.)
        if espressomd.has_features(["VIRTUAL_SITES_INERTIALESS_TRACERS"]):
            volcons = espressomd.interactions.IBM_VolCons(softID=1, kappaV=1.)
        system = self.system
        system.bonded_inter.add(angle)
        system.bonded_inter.add(dihe)
        if espressomd.has_features(["VIRTUAL_SITES_INERTIALESS_TRACERS"]):
            system.bonded_inter.add(volcons)

        self.system.part.add(pos=[0.1, 0.1, 0.1])
        self.system.part.add(pos=[0.2, 0.1, 0.1])
        self.system.part.add(pos=[0.2, 0.2, 0.1])
        self.system.part.add(pos=[0.1, 0.2, 0.1])
        system.part.by_id(0).add_bond((angle, 1, 2))
        system.part.by_id(0).add_bond((dihe, 1, 2, 3))
        if espressomd.has_features(["VIRTUAL_SITES_INERTIALESS_TRACERS"]):
            system.part.by_id(0).add_bond((volcons,))

        system.auto_exclusions(distance=2)

        for p in system.part.all():
            self.assertEqual(list(p.exclusions), [])


if __name__ == "__main__":
    ut.main()
