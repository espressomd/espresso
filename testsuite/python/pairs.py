#
# Copyright (C) 2010-2022 The ESPResSo project
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


class PairTest(ut.TestCase):
    """
    Tests the particle pair finder of the cell system.
    It checks that particles are found if their distance is below the threshold,
    no matter the type of cell system, periodicity and the image box they are in.
    Also tests that the ``types`` argument works as expected and an exception is raised
    when the distance threshold is larger than the cell size.
    """

    system = espressomd.System(box_l=3 * [10.])
    system.time_step = 0.1
    system.cell_system.skin = 0.3
    # Force an appropriate cell grid
    system.min_global_cut = 1.6

    def setUp(self):
        vel = [1., 2., 3.]

        self.system.part.add(id=0, pos=[5., 4.5, 5.], v=vel, type=0)
        self.system.part.add(id=1, pos=[5., 5.5, 5.], v=vel, type=1)
        self.system.part.add(id=2, pos=[9.5, 5., 5.], v=vel, type=2)
        self.system.part.add(id=3, pos=[0.5, 5., 5.], v=vel, type=3)
        self.system.part.add(id=4, pos=[5., 5., 9.5], v=vel, type=4)
        self.system.part.add(id=5, pos=[5., 5., 0.5], v=vel, type=5)
        self.system.part.add(id=6, pos=[5., 9.5, 5.], v=vel, type=6)
        self.system.part.add(id=7, pos=[5., 0.5, 5.], v=vel, type=7)
        self.system.part.add(id=8, pos=[5., 9.5, 9.5], v=vel, type=8)
        self.system.part.add(id=9, pos=[5., 0.5, 0.5], v=vel, type=9)
        self.system.part.add(id=10, pos=[1., 1., 1.], v=vel, type=10)
        self.system.part.add(id=11, pos=[9., 9., 9.], v=vel, type=11)

        self.types_to_get_pairs = [0, 1, 4, 5, 6]

    def tearDown(self):
        self.system.part.clear()

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
        pairs = self.system.cell_system.get_pairs(1.5)
        epairs = self.expected_pairs(self.system.periodicity)
        self.assertSetEqual(set(pairs), set(epairs))

        pairs_by_type = self.system.cell_system.get_pairs(
            1.5, types=self.types_to_get_pairs)
        epairs_by_type = self.expected_pairs_with_types(
            self.system.periodicity)
        self.assertSetEqual(set(pairs_by_type), set(epairs_by_type))

    def test_input_exceptions(self):
        with self.assertRaisesRegex(ValueError, "Unknown argument types='none'"):
            self.system.cell_system.get_pairs(0.1, types="none")
        with self.assertRaisesRegex(RuntimeError, "Provided argument of type 'int' is not convertible to 'std::vector<int>'"):
            self.system.cell_system.get_pairs(0.1, types=3)
        with self.assertRaisesRegex(RuntimeError, "Provided argument of type .+ because it contains a value that is not convertible to 'int'"):
            self.system.cell_system.get_pairs(0.1, types={'3.': 6.})
        # check no exception for list of length 1
        self.system.cell_system.get_pairs(0.1, types=[3])

    def check_range_exception(self):
        with self.assertRaisesRegex(ValueError, "pair search distance 3.* bigger than the decomposition range"):
            self.system.cell_system.get_pairs(3.)

    def run_and_check(self, n_steps=100):
        self.system.integrator.run(0)
        self.check_pairs()
        self.system.integrator.run(n_steps)
        self.check_pairs()

    def test_nsquare(self):
        self.system.cell_system.set_n_square()
        self.system.periodicity = [True, True, True]
        self.run_and_check()

    def test_nsquare_partial_z(self):
        self.system.cell_system.set_n_square()
        self.system.periodicity = [True, True, False]
        self.run_and_check()

    def test_dd(self):
        self.system.cell_system.set_regular_decomposition()
        self.system.periodicity = [True, True, True]
        self.run_and_check()
        self.check_range_exception()

    def test_dd_partial_z(self):
        self.system.cell_system.set_regular_decomposition()
        self.system.periodicity = [True, True, False]
        self.run_and_check()
        self.check_range_exception()

    def test_non_consecutive(self):
        self.system.part.clear()
        self.system.part.add(id=100, pos=(0.1, 0.5, 0.5))
        self.system.part.add(id=200, pos=(0.2, 0.5, 0.5))
        self.system.part.add(id=1, pos=(0.3, 0.5, 0.5))
        self.system.part.add(id=2, pos=(0.4, 0.5, 0.5))

        self.system.integrator.run(0)

        pairs = self.system.part.pairs()
        pairs_ids = []
        for pair in pairs:
            pairs_ids.append(tuple(sorted([pair[0].id, pair[1].id])))
        expected_pairs = [(1, 2), (1, 100), (2, 200),
                          (100, 200), (2, 100), (1, 200)]
        self.assertSetEqual(set(pairs_ids), set(expected_pairs))


if __name__ == "__main__":
    ut.main()
