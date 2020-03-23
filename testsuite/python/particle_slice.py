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
import unittest_decorators as utx
import espressomd
from espressomd import has_features
import numpy as np


class ParticleSliceTest(ut.TestCase):

    state = [[0, 0, 0], [0, 0, 1]]
    system = espressomd.System(box_l=[10, 10, 10])

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.system.part.clear()
        for i in range(4):
            self.system.part.add(pos=[0, 0, i])

        if has_features(["EXTERNAL_FORCES"]):
            self.system.part[1].fix = self.state[1]
            np.testing.assert_array_equal(
                np.copy(self.system.part[0].fix), self.state[0])
            np.testing.assert_array_equal(
                np.copy(self.system.part[1].fix), self.state[1])
            np.testing.assert_array_equal(
                np.copy(self.system.part[:2].fix), self.state)
        xs = self.system.part[:].pos
        for i in range(len(xs)):
            np.testing.assert_array_equal(
                xs[i], np.copy(self.system.part[i].pos))

    @utx.skipIfMissingFeatures(["EXTERNAL_FORCES"])
    def test_1_set_different_values(self):
        self.state[0] = [1, 0, 0]
        self.state[1] = [1, 0, 0]
        self.system.part[:2].fix = self.state
        np.testing.assert_array_equal(
            np.copy(self.system.part[:2].fix), self.state)

    @utx.skipIfMissingFeatures(["EXTERNAL_FORCES"])
    def test_2_set_same_value(self):
        self.state[0] = [0, 1, 0]
        self.state[1] = [0, 1, 0]
        self.system.part[:2].fix = self.state[1]
        np.testing.assert_array_equal(
            np.copy(self.system.part[:2].fix), self.state)

    @utx.skipIfMissingFeatures(["EXTERNAL_FORCES"])
    def test_3_set_one_value(self):
        self.state[1] = [0, 0, 1]
        self.system.part[1:2].fix = self.state[1]
        np.testing.assert_array_equal(
            np.copy(self.system.part[:2].fix), self.state)

    @utx.skipIfMissingFeatures(["EXTERNAL_FORCES"])
    def test_4_str(self):
        self.assertEqual(repr(self.system.part[0].fix),
                         repr(np.array([0, 1, 0])))
        self.assertEqual(repr(self.system.part[:2].fix),
                         repr(np.array([[0, 1, 0], [0, 0, 1]])))

    def test_pos_str(self):
        self.system.part[0].pos = [0, 0, 0]
        self.system.part[1].pos = [0, 0, 1]
        self.assertEqual(repr(self.system.part[0].pos),
                         repr(np.array([0.0, 0.0, 0.0])))
        self.assertEqual(repr(self.system.part[:2].pos),
                         repr(np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])))

    @utx.skipIfMissingFeatures(["ELECTROSTATICS"])
    def test_scalar(self):
        self.system.part[:1].q = 1.3
        self.assertEqual(self.system.part[0].q, 1.3)
        self.system.part[:2].q = 2.0
        self.assertEqual(self.system.part[0].q, 2)
        self.assertEqual(self.system.part[1].q, 2)
        self.system.part[:2].q = 3
        self.assertEqual(self.system.part[0].q, 3)
        self.assertEqual(self.system.part[1].q, 3)
        self.system.part[:2].q = [-1, 1.0]
        self.assertEqual(self.system.part[0].q, -1)
        self.assertEqual(self.system.part[1].q, 1)
        qs = self.system.part[:2].q
        self.assertEqual(qs[0], -1)
        self.assertEqual(qs[1], 1)

    def test_bonds(self):

        fene = espressomd.interactions.FeneBond(k=1, d_r_max=1, r_0=1)
        self.system.bonded_inter.add(fene)

        # Setter

        # tuple
        b = fene, 0
        self.system.part[2].bonds = b
        self.assertEqual(self.system.part[:].bonds, [(), (), ((fene, 0),), ()])

        # list
        self.system.part[:].bonds = []
        b = [fene, 0]
        self.system.part[2].bonds = b
        self.assertEqual(self.system.part[:].bonds, [(), (), ((fene, 0),), ()])

        # nested list single
        self.system.part[:].bonds = []
        b = [[fene, 0]]
        self.system.part[2].bonds = b
        self.assertEqual(self.system.part[:].bonds, [(), (), ((fene, 0),), ()])

        # nested list multi
        self.system.part[:].bonds = []
        b = [[fene, 0], [fene, 1]]
        self.system.part[2].bonds = b
        self.assertEqual(self.system.part[:].bonds,
                         [(), (), ((fene, 0), (fene, 1)), ()])

        # nested tuple single
        self.system.part[:].bonds = []
        b = ((fene, 0),)
        self.system.part[2].bonds = b
        self.assertEqual(self.system.part[:].bonds, [(), (), ((fene, 0),), ()])

        # nested tuple multi
        self.system.part[:].bonds = []
        b = ((fene, 0), (fene, 1))
        self.system.part[2].bonds = b
        self.assertEqual(self.system.part[:].bonds,
                         [(), (), ((fene, 0), (fene, 1)), ()])

        # Add/Del bonds
        self.system.part[:].bonds = []
        self.system.part[2].add_bond((fene, 0))
        self.assertEqual(self.system.part[:].bonds, [(), (), ((fene, 0),), ()])
        self.system.part[2].add_bond((fene, 1))
        self.assertEqual(self.system.part[:].bonds,
                         [(), (), ((fene, 0), (fene, 1)), ()])
        self.system.part[2].delete_bond((fene, 1))
        self.assertEqual(self.system.part[:].bonds, [(), (), ((fene, 0),), ()])
        self.system.part[2].delete_bond((fene, 0))
        self.assertEqual(self.system.part[:].bonds, [(), (), (), ()])

        self.system.part[:].bonds = []
        self.system.part[2].add_bond((fene, 0))
        self.system.part[2].add_bond((fene, 1))
        self.system.part[2].delete_all_bonds()
        self.assertEqual(self.system.part[:].bonds, [(), (), (), ()])

        # Slices

        # tuple for all
        self.system.part[:].bonds = []
        b = fene, 0
        self.system.part[2:].bonds = b
        self.assertEqual(self.system.part[:].bonds,
                         [(), (), ((fene, 0),), ((fene, 0),)])

        # list for all
        self.system.part[:].bonds = []
        b = [fene, 0]
        self.system.part[2:].bonds = b
        self.assertEqual(self.system.part[:].bonds,
                         [(), (), ((fene, 0),), ((fene, 0),)])

        # nested list single for all
        self.system.part[:].bonds = []
        b = [[fene, 0]]
        self.system.part[2:].bonds = b
        self.assertEqual(self.system.part[:].bonds,
                         [(), (), ((fene, 0),), ((fene, 0),)])

        # nested list multi for all
        self.system.part[:].bonds = []
        b = [[fene, 0], [fene, 1]]
        self.system.part[2:].bonds = b
        self.assertEqual(self.system.part[:].bonds,
                         [(), (), ((fene, 0), (fene, 1)), ((fene, 0), (fene, 1))])

        # tuples for each
        self.system.part[:].bonds = []
        b = (((fene, 0),), ((fene, 1),))
        self.system.part[2:].bonds = b
        self.assertEqual(self.system.part[:].bonds,
                         [(), (), ((fene, 0),), ((fene, 1),)])

        # lists for each
        self.system.part[:].bonds = []
        b = [[[fene, 0]], [[fene, 1]]]
        self.system.part[2:].bonds = b
        self.assertEqual(self.system.part[:].bonds,
                         [(), (), ((fene, 0),), ((fene, 1),)])

        # multi tuples for each
        self.system.part[:].bonds = []
        b = (((fene, 0), (fene, 1)), ((fene, 0), (fene, 1)))
        self.system.part[2:].bonds = b
        self.assertEqual(self.system.part[:].bonds,
                         [(), (), ((fene, 0), (fene, 1)), ((fene, 0), (fene, 1))])

        # multi lists for each
        self.system.part[:].bonds = []
        b = [[[fene, 0], [fene, 1]], [[fene, 0], [fene, 1]]]
        self.system.part[2:].bonds = b
        self.assertEqual(self.system.part[:].bonds,
                         [(), (), ((fene, 0), (fene, 1)), ((fene, 0), (fene, 1))])

        # Add/Del bonds
        self.system.part[:].bonds = []
        self.system.part[2:].add_bond((fene, 0))
        self.assertEqual(self.system.part[:].bonds,
                         [(), (), ((fene, 0),), ((fene, 0),)])
        self.system.part[2:].add_bond((fene, 1))
        self.assertEqual(self.system.part[:].bonds,
                         [(), (), ((fene, 0), (fene, 1)), ((fene, 0), (fene, 1))])
        self.system.part[2:].delete_bond((fene, 1))
        self.assertEqual(self.system.part[:].bonds,
                         [(), (), ((fene, 0),), ((fene, 0),)])
        self.system.part[2:].delete_bond((fene, 0))
        self.assertEqual(self.system.part[:].bonds, [(), (), (), ()])

        b = [[[fene, 0], [fene, 1]], [[fene, 0], [fene, 1]]]
        self.system.part[2:].bonds = b
        self.system.part[:].delete_all_bonds()
        self.assertEqual(self.system.part[:].bonds, [(), (), (), ()])

    @utx.skipIfMissingFeatures(["EXCLUSIONS"])
    def test_exclusions(self):

        def assert_exclusions_equal(B):
            A = self.system.part[:].exclusions
            self.assertEqual([x.tolist() for x in A], B)

        # Setter
        # int
        self.system.part[:].exclusions = []
        b = 1
        self.system.part[2].exclusions = b
        assert_exclusions_equal([[], [2], [1], []])

        # single list
        self.system.part[:].exclusions = []
        b = [1]
        self.system.part[2].exclusions = b
        assert_exclusions_equal([[], [2], [1], []])

        # tuple
        self.system.part[:].exclusions = []
        b = (0, 1)
        self.system.part[2].exclusions = b
        assert_exclusions_equal([[2], [2], [0, 1], []])

        # list
        self.system.part[:].exclusions = []
        b = [0, 1]
        self.system.part[2].exclusions = b
        assert_exclusions_equal([[2], [2], [0, 1], []])

        # Add/Del exclusions
        self.system.part[:].exclusions = []
        self.system.part[2].add_exclusion(1)
        assert_exclusions_equal([[], [2], [1], []])
        self.system.part[2].add_exclusion(0)
        assert_exclusions_equal([[2], [2], [1, 0], []])
        self.system.part[2].delete_exclusion(0)
        assert_exclusions_equal([[], [2], [1], []])
        self.system.part[2].delete_exclusion(1)
        assert_exclusions_equal([[], [], [], []])

        # Slices

        # single list for all
        self.system.part[:].exclusions = []
        b = [1]
        self.system.part[2:].exclusions = b
        assert_exclusions_equal([[], [2, 3], [1], [1]])

        # list for all
        self.system.part[:].exclusions = []
        b = [0, 1]
        self.system.part[2:].exclusions = b
        assert_exclusions_equal([[2, 3], [2, 3], [0, 1], [0, 1]])

        # single list for each
        self.system.part[:].exclusions = []
        b = [[0], [0]]
        self.system.part[2:].exclusions = b
        assert_exclusions_equal([[2, 3], [], [0], [0]])

        # multi list for each
        self.system.part[:].exclusions = []
        b = [[0, 1], [0, 1]]
        self.system.part[2:].exclusions = b
        assert_exclusions_equal([[2, 3], [2, 3], [0, 1], [0, 1]])

        # Add/Del exclusions
        self.system.part[:].exclusions = []
        self.system.part[2:].add_exclusion(1)
        assert_exclusions_equal([[], [2, 3], [1], [1]])
        self.system.part[2:].add_exclusion(0)
        assert_exclusions_equal([[2, 3], [2, 3], [1, 0], [1, 0]])
        self.system.part[2:].delete_exclusion(0)
        assert_exclusions_equal([[], [2, 3], [1], [1]])
        self.system.part[2:].delete_exclusion(1)
        assert_exclusions_equal([[], [], [], []])

    @utx.skipIfMissingFeatures("VIRTUAL_SITES_RELATIVE")
    def test_vs_relative(self):

        self.system.part.clear()
        self.system.part.add(pos=[0, 0, 0])
        self.system.part.add(pos=[0, 0, 0])
        self.system.part.add(pos=[0, 0, 0])
        self.system.part.add(pos=[0, 0, 0])
        self.system.part[0].vs_relative = [1, 1.0, (1.0, 1.0, 1.0, 1.0)]

        self.assertEqual(repr(self.system.part[:].vs_relative),
                         repr([(1, 1.0, np.array([1., 1., 1., 1.])),
                               (0, 0.0, np.array([0., 0., 0., 0.])),
                               (0, 0.0, np.array([0., 0., 0., 0.])),
                               (0, 0.0, np.array([0., 0., 0., 0.]))]))

        self.system.part[:].vs_relative = [1, 1.0, (1.0, 1.0, 1.0, 1.0)]

        self.assertEqual(repr(self.system.part[:].vs_relative),
                         repr([(1, 1.0, np.array([1., 1., 1., 1.])),
                               (1, 1.0, np.array([1., 1., 1., 1.])),
                               (1, 1.0, np.array([1., 1., 1., 1.])),
                               (1, 1.0, np.array([1., 1., 1., 1.]))]))

        self.system.part[:].vs_relative = [[1, 1.0, (1.0, 1.0, 1.0, 1.0)],
                                           [1, 2.0, (1.0, 1.0, 1.0, 1.0)],
                                           [1, 3.0, (1.0, 1.0, 1.0, 1.0)],
                                           [1, 4.0, (1.0, 1.0, 1.0, 1.0)]]

        self.assertEqual(repr(self.system.part[:].vs_relative),
                         repr([(1, 1.0, np.array([1., 1., 1., 1.])),
                               (1, 2.0, np.array([1., 1., 1., 1.])),
                               (1, 3.0, np.array([1., 1., 1., 1.])),
                               (1, 4.0, np.array([1., 1., 1., 1.]))]))

    def test_multiadd(self):
        self.system.part.clear()
        positions = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        self.system.part.add(pos=positions)
        for p in positions:
            self.system.part.add(pos=p)
        np.testing.assert_allclose(
            np.copy(self.system.part[:3].pos), np.copy(self.system.part[3:6].pos))

        with self.assertRaises(ValueError):
            self.system.part.add(pos=([1, 1, 1], [2, 2, 2]), type=0)
        with self.assertRaises(ValueError):
            self.system.part.add(pos=([1, 1, 1], [2, 2, 2]), type=(0, 1, 2))

        self.system.part.clear()
        self.system.part.add(pos=([1, 1, 1], [2, 2, 2]), type=(0, 1))
        self.assertEqual(self.system.part[0].type, 0)
        self.assertEqual(self.system.part[1].type, 1)

    def test_empty(self):
        np.testing.assert_array_equal(self.system.part[0:0].pos, np.empty(0))

    def test_len(self):
        self.assertEqual(len(self.system.part[0:0]), 0)
        self.assertEqual(len(self.system.part[0:1]), 1)
        self.assertEqual(len(self.system.part[0:2]), 2)

    def test_non_existing_property(self):
        with self.assertRaises(AttributeError):
            self.system.part[:].thispropertydoesnotexist = 1.0


if __name__ == "__main__":
    ut.main()
