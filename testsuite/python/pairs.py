# Copyright (C) 2010-2019 The ESPResSo project
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
import espressomd


class PairTest(ut.TestCase):
    """
    Tests the particle pair finder of the cell system.
    It checks that particles are found if their distance is below the threshold, 
    no matter the type of cell system, periodicity and the image box they are in.
    Also tests that the ``types`` argument works as expected and an exception is raised
    when the distance threshold is larger than the cell size.
    """

    s = espressomd.System(box_l=[1.0, 1.0, 1.0])

    def setUp(self):
        self.s.time_step = 0.1     
        self.s.box_l = 3 * [10.]
        self.s.cell_system.skin = 0.3

        # Force an appropriate cell grid
        self.s.min_global_cut = 1.6

        vel = [1., 2., 3.]

        self.s.part.add(id=0, pos=[5., 4.5, 5.], v=vel, type=0)
        self.s.part.add(id=1, pos=[5., 5.5, 5.], v=vel, type=1)
        self.s.part.add(id=2, pos=[9.5, 5., 5.], v=vel, type=2)
        self.s.part.add(id=3, pos=[0.5, 5., 5.], v=vel, type=3)
        self.s.part.add(id=4, pos=[5., 5., 9.5], v=vel, type=4)
        self.s.part.add(id=5, pos=[5., 5., 0.5], v=vel, type=5)
        self.s.part.add(id=6, pos=[5., 9.5, 5.], v=vel, type=6)
        self.s.part.add(id=7, pos=[5., 0.5, 5.], v=vel, type=7)
        self.s.part.add(id=8, pos=[5., 9.5, 9.5], v=vel, type=8)
        self.s.part.add(id=9, pos=[5., 0.5, 0.5], v=vel, type=9)
        self.s.part.add(id=10, pos=[1., 1., 1.], v=vel, type=10)
        self.s.part.add(id=11, pos=[9., 9., 9.], v=vel, type=11)

        self.types_to_get_pairs = [0, 1, 4, 5, 6]

    def tearDown(self):
        self.s.part.clear()

    def expected_pairs(self, periodicity):
        if all(periodicity == (1, 1, 1)):
            return [(0, 1), (2, 3), (4, 5), (6, 7), (8, 9)]
        elif all(periodicity == (1, 1, 0)):
            return [(0, 1), (2, 3), (6, 7)]

    def expected_pairs_with_types(self, periodicity):
        if all(periodicity == (1, 1, 1)):
            return [(0, 1), (4, 5)]
        elif all(periodicity == (1, 1, 0)):
            return [(0, 1)]

    def check_pairs(self):
        pairs = self.s.cell_system.get_pairs(1.5)
        epairs = self.expected_pairs(self.s.periodicity)
        self.assertSetEqual(set(pairs), set(epairs))

        pairs_by_type = self.s.cell_system.get_pairs(
            1.5, types=self.types_to_get_pairs)
        epairs_by_type = self.expected_pairs_with_types(self.s.periodicity)
        self.assertSetEqual(set(pairs_by_type), set(epairs_by_type))

    def test_input_exceptions(self):
        with self.assertRaises(ValueError):
            self.s.cell_system.get_pairs(0.1, types=3)
        # check no exception for list of length 1
        self.s.cell_system.get_pairs(0.1, types=[3])

    def check_range_exception(self):
        with self.assertRaises(Exception):
            self.s.cell_system.get_pairs(3.)

    def run_and_check(self, n_steps=100):
        self.s.integrator.run(0)
        self.check_pairs()
        self.s.integrator.run(n_steps)
        self.check_pairs()

    def test_nsquare(self):
        self.s.cell_system.set_n_square()
        self.s.periodicity = [1, 1, 1]
        self.run_and_check()

    def test_nsquare_partial_z(self):
        self.s.cell_system.set_n_square()
        self.s.periodicity = [1, 1, 0]
        self.run_and_check()

    def test_dd(self):
        self.s.cell_system.set_domain_decomposition()
        self.s.periodicity = [1, 1, 1]
        self.run_and_check()
        self.check_range_exception()

    def test_dd_partial_z(self):
        self.s.cell_system.set_domain_decomposition()
        self.s.periodicity = [1, 1, 0]
        self.run_and_check()
        self.check_range_exception()


if __name__ == "__main__":
    ut.main()
