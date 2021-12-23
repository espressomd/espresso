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
import numpy as np


class ParticleSliceTest(ut.TestCase):

    system = espressomd.System(box_l=3 * [10])

    def setUp(self):
        self.fix_flags = [[False, False, False], [False, False, True]]
        self.p0 = self.system.part.add(pos=[0, 0, 0])
        self.p1 = self.system.part.add(pos=[0, 0, 1])
        self.p2 = self.system.part.add(pos=[0, 0, 2])
        self.p3 = self.system.part.add(pos=[0, 0, 3])

        self.all_partcls = self.system.part.all()
        self.p0p1 = self.system.part.by_ids([0, 1])
        self.p2p3 = self.system.part.by_ids([2, 3])

        if espressomd.has_features(["EXTERNAL_FORCES"]):
            self.p1.fix = self.fix_flags[1]
            np.testing.assert_array_equal(
                np.copy(self.p0.fix), self.fix_flags[0])
            np.testing.assert_array_equal(
                np.copy(self.p1.fix), self.fix_flags[1])
            np.testing.assert_array_equal(
                np.copy(self.p0p1.fix), self.fix_flags)

    def tearDown(self):
        self.system.part.clear()

    @utx.skipIfMissingFeatures(["EXTERNAL_FORCES"])
    def test_1_set_different_values(self):
        self.fix_flags[0] = [1, 0, 0]
        self.fix_flags[1] = [1, 0, 0]
        self.p0p1.fix = self.fix_flags
        np.testing.assert_array_equal(
            np.copy(self.p0p1.fix), self.fix_flags)

    @utx.skipIfMissingFeatures(["EXTERNAL_FORCES"])
    def test_2_set_same_value(self):
        self.fix_flags[0] = [0, 1, 0]
        self.fix_flags[1] = [0, 1, 0]
        self.p0p1.fix = self.fix_flags[1]
        np.testing.assert_array_equal(
            np.copy(self.p0p1.fix), self.fix_flags)

    @utx.skipIfMissingFeatures(["EXTERNAL_FORCES"])
    def test_3_set_one_value(self):
        self.fix_flags[1] = [0, 0, 1]
        self.p1.fix = self.fix_flags[1]
        np.testing.assert_array_equal(
            np.copy(self.p0p1.fix), self.fix_flags)

    @utx.skipIfMissingFeatures(["EXTERNAL_FORCES"])
    def test_4_str(self):
        self.assertEqual(repr(self.p0.fix),
                         repr(np.array([0, 0, 0])))
        self.assertEqual(repr(self.p0p1.fix),
                         repr(np.array([[0, 0, 0], [0, 0, 1]])))

    def test_pos_str(self):
        self.p0.pos = [0, 0, 0]
        self.p1.pos = [0, 0, 1]
        self.assertEqual(repr(self.p0.pos),
                         repr(np.array([0.0, 0.0, 0.0])))
        self.assertEqual(repr(self.p0p1.pos),
                         repr(np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])))

    @utx.skipIfMissingFeatures(["ELECTROSTATICS"])
    def test_scalar(self):
        self.p0.q = 1.3
        self.assertEqual(self.p0.q, 1.3)
        self.p0p1.q = 2.0
        self.assertEqual(self.p0.q, 2)
        self.assertEqual(self.p1.q, 2)
        self.p0p1.q = 3
        self.assertEqual(self.p0.q, 3)
        self.assertEqual(self.p1.q, 3)
        self.p0p1.q = [-1, 1.0]
        self.assertEqual(self.p0.q, -1)
        self.assertEqual(self.p1.q, 1)
        qs = self.p0p1.q
        self.assertEqual(qs[0], -1)
        self.assertEqual(qs[1], 1)

    def test_bonds(self):

        fene = espressomd.interactions.FeneBond(k=1, d_r_max=1, r_0=1)
        self.system.bonded_inter.add(fene)

        # Setter

        # tuple
        b = fene, 0
        self.p2.bonds = b
        self.assertEqual(self.all_partcls.bonds, [(), (), ((fene, 0),), ()])

        # list
        self.all_partcls.bonds = []
        b = [fene, 0]
        self.p2.bonds = b
        self.assertEqual(self.all_partcls.bonds, [(), (), ((fene, 0),), ()])

        # nested list single
        self.all_partcls.bonds = []
        b = [[fene, 0]]
        self.p2.bonds = b
        self.assertEqual(self.all_partcls.bonds, [(), (), ((fene, 0),), ()])

        # nested list multi
        self.all_partcls.bonds = []
        b = [[fene, 0], [fene, 1]]
        self.p2.bonds = b
        self.assertEqual(self.all_partcls.bonds,
                         [(), (), ((fene, 0), (fene, 1)), ()])

        # nested tuple single
        self.all_partcls.bonds = []
        b = ((fene, 0),)
        self.p2.bonds = b
        self.assertEqual(self.all_partcls.bonds, [(), (), ((fene, 0),), ()])

        # nested tuple multi
        self.all_partcls.bonds = []
        b = ((fene, 0), (fene, 1))
        self.p2.bonds = b
        self.assertEqual(self.all_partcls.bonds,
                         [(), (), ((fene, 0), (fene, 1)), ()])

        # Add/Del bonds
        self.all_partcls.bonds = []
        self.p2.add_bond((fene, 0))
        self.assertEqual(self.all_partcls.bonds, [(), (), ((fene, 0),), ()])
        self.p2.add_bond((fene, 1))
        self.assertEqual(self.all_partcls.bonds,
                         [(), (), ((fene, 0), (fene, 1)), ()])
        self.p2.delete_bond((fene, 1))
        self.assertEqual(self.all_partcls.bonds, [(), (), ((fene, 0),), ()])
        self.p2.delete_bond((fene, 0))
        self.assertEqual(self.all_partcls.bonds, [(), (), (), ()])

        self.all_partcls.bonds = []
        self.p2.add_bond((fene, 0))
        self.p2.add_bond((fene, 1))
        self.p2.delete_all_bonds()
        self.assertEqual(self.all_partcls.bonds, [(), (), (), ()])

        # Slices

        # tuple for all
        self.all_partcls.bonds = []
        b = fene, 0
        self.p2p3.bonds = b
        self.assertEqual(self.all_partcls.bonds,
                         [(), (), ((fene, 0),), ((fene, 0),)])

        # list for all
        self.all_partcls.bonds = []
        b = [fene, 0]
        self.p2p3.bonds = b
        self.assertEqual(self.all_partcls.bonds,
                         [(), (), ((fene, 0),), ((fene, 0),)])

        # nested list single for all
        self.all_partcls.bonds = []
        b = [[fene, 0]]
        self.p2p3.bonds = b
        self.assertEqual(self.all_partcls.bonds,
                         [(), (), ((fene, 0),), ((fene, 0),)])

        # nested list multi for all
        self.all_partcls.bonds = []
        b = [[fene, 0], [fene, 1]]
        self.p2p3.bonds = b
        self.assertEqual(self.all_partcls.bonds,
                         [(), (), ((fene, 0), (fene, 1)), ((fene, 0), (fene, 1))])

        # tuples for each
        self.all_partcls.bonds = []
        b = (((fene, 0),), ((fene, 1),))
        self.p2p3.bonds = b
        self.assertEqual(self.all_partcls.bonds,
                         [(), (), ((fene, 0),), ((fene, 1),)])

        # lists for each
        self.all_partcls.bonds = []
        b = [[[fene, 0]], [[fene, 1]]]
        self.p2p3.bonds = b
        self.assertEqual(self.all_partcls.bonds,
                         [(), (), ((fene, 0),), ((fene, 1),)])

        # multi tuples for each
        self.all_partcls.bonds = []
        b = (((fene, 0), (fene, 1)), ((fene, 0), (fene, 1)))
        self.p2p3.bonds = b
        self.assertEqual(self.all_partcls.bonds,
                         [(), (), ((fene, 0), (fene, 1)), ((fene, 0), (fene, 1))])

        # multi lists for each
        self.all_partcls.bonds = []
        b = [[[fene, 0], [fene, 1]], [[fene, 0], [fene, 1]]]
        self.p2p3.bonds = b
        self.assertEqual(self.all_partcls.bonds,
                         [(), (), ((fene, 0), (fene, 1)), ((fene, 0), (fene, 1))])

        # Add/Del bonds
        self.all_partcls.bonds = []
        self.p2p3.add_bond((fene, 0))
        self.assertEqual(self.all_partcls.bonds,
                         [(), (), ((fene, 0),), ((fene, 0),)])
        self.p2p3.add_bond((fene, 1))
        self.assertEqual(self.all_partcls.bonds,
                         [(), (), ((fene, 0), (fene, 1)), ((fene, 0), (fene, 1))])
        self.p2p3.delete_bond((fene, 1))
        self.assertEqual(self.all_partcls.bonds,
                         [(), (), ((fene, 0),), ((fene, 0),)])
        self.p2p3.delete_bond((fene, 0))
        self.assertEqual(self.all_partcls.bonds, [(), (), (), ()])

        b = [[[fene, 0], [fene, 1]], [[fene, 0], [fene, 1]]]
        self.p2p3.bonds = b
        self.all_partcls.delete_all_bonds()
        self.assertEqual(self.all_partcls.bonds, [(), (), (), ()])

    @utx.skipIfMissingFeatures(["EXCLUSIONS"])
    def test_exclusions(self):

        def assert_exclusions_equal(B):
            A = self.all_partcls.exclusions
            self.assertEqual([x.tolist() for x in A], B)

        # Setter
        # int
        self.all_partcls.exclusions = []
        b = 1
        self.p2.exclusions = b
        assert_exclusions_equal([[], [2], [1], []])

        # single list
        self.all_partcls.exclusions = []
        b = [1]
        self.p2.exclusions = b
        assert_exclusions_equal([[], [2], [1], []])

        # tuple
        self.all_partcls.exclusions = []
        b = (0, 1)
        self.p2.exclusions = b
        assert_exclusions_equal([[2], [2], [0, 1], []])

        # list
        self.all_partcls.exclusions = []
        b = [0, 1]
        self.p2.exclusions = b
        assert_exclusions_equal([[2], [2], [0, 1], []])

        # Add/Del exclusions
        self.all_partcls.exclusions = []
        self.p2.add_exclusion(1)
        assert_exclusions_equal([[], [2], [1], []])
        self.p2.add_exclusion(0)
        assert_exclusions_equal([[2], [2], [1, 0], []])
        self.p2.delete_exclusion(0)
        assert_exclusions_equal([[], [2], [1], []])
        self.p2.delete_exclusion(1)
        assert_exclusions_equal([[], [], [], []])

        # Slices

        # single list for all
        self.all_partcls.exclusions = []
        b = [1]
        self.p2p3.exclusions = b
        assert_exclusions_equal([[], [2, 3], [1], [1]])

        # list for all
        self.all_partcls.exclusions = []
        b = [0, 1]
        self.p2p3.exclusions = b
        assert_exclusions_equal([[2, 3], [2, 3], [0, 1], [0, 1]])

        # single list for each
        self.all_partcls.exclusions = []
        b = [[0], [0]]
        self.p2p3.exclusions = b
        assert_exclusions_equal([[2, 3], [], [0], [0]])

        # multi list for each
        self.all_partcls.exclusions = []
        b = [[0, 1], [0, 1]]
        self.p2p3.exclusions = b
        assert_exclusions_equal([[2, 3], [2, 3], [0, 1], [0, 1]])

        # Add/Del exclusions
        self.all_partcls.exclusions = []
        self.p2p3.add_exclusion(1)
        assert_exclusions_equal([[], [2, 3], [1], [1]])
        self.p2p3.add_exclusion(0)
        assert_exclusions_equal([[2, 3], [2, 3], [1, 0], [1, 0]])
        self.p2p3.delete_exclusion(0)
        assert_exclusions_equal([[], [2, 3], [1], [1]])
        self.p2p3.delete_exclusion(1)
        assert_exclusions_equal([[], [], [], []])

    @utx.skipIfMissingFeatures("VIRTUAL_SITES_RELATIVE")
    def test_vs_relative(self):

        self.system.part.clear()
        p0 = self.system.part.add(pos=[0, 0, 0])
        self.system.part.add(pos=3 * [[0, 0, 0]])
        p0.vs_relative = [1, 1.0, (1.0, 1.0, 1.0, 1.0)]
        all_partcls = self.system.part.all()

        self.assertEqual(repr(all_partcls.vs_relative),
                         repr([(1, 1.0, np.array([1., 1., 1., 1.])),
                               (0, 0.0, np.array([1., 0., 0., 0.])),
                               (0, 0.0, np.array([1., 0., 0., 0.])),
                               (0, 0.0, np.array([1., 0., 0., 0.]))]))

        all_partcls.vs_relative = [1, 1.0, (1.0, 1.0, 1.0, 1.0)]

        self.assertEqual(repr(all_partcls.vs_relative),
                         repr([(1, 1.0, np.array([1., 1., 1., 1.])),
                               (1, 1.0, np.array([1., 1., 1., 1.])),
                               (1, 1.0, np.array([1., 1., 1., 1.])),
                               (1, 1.0, np.array([1., 1., 1., 1.]))]))

        all_partcls.vs_relative = [[1, 1.0, (1.0, 1.0, 1.0, 1.0)],
                                   [1, 2.0, (1.0, 1.0, 1.0, 1.0)],
                                   [1, 3.0, (1.0, 1.0, 1.0, 1.0)],
                                   [1, 4.0, (1.0, 1.0, 1.0, 1.0)]]

        self.assertEqual(repr(all_partcls.vs_relative),
                         repr([(1, 1.0, np.array([1., 1., 1., 1.])),
                               (1, 2.0, np.array([1., 1., 1., 1.])),
                               (1, 3.0, np.array([1., 1., 1., 1.])),
                               (1, 4.0, np.array([1., 1., 1., 1.]))]))

    def test_multiadd(self):
        self.system.part.clear()
        positions = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        partcls_multiadd = self.system.part.add(pos=positions)
        for p in positions:
            self.system.part.add(pos=p)
        partcls_singleadd = self.system.part.by_ids(range(3, 6))
        np.testing.assert_allclose(
            np.copy(partcls_multiadd.pos), np.copy(partcls_singleadd.pos))

        with self.assertRaises(ValueError):
            self.system.part.add(pos=([1, 1, 1], [2, 2, 2]), type=0)
        with self.assertRaises(ValueError):
            self.system.part.add(pos=([1, 1, 1], [2, 2, 2]), type=(0, 1, 2))

        self.system.part.clear()
        p0, p1 = self.system.part.add(pos=([1, 1, 1], [2, 2, 2]), type=(0, 1))
        self.assertEqual(p0.type, 0)
        self.assertEqual(p1.type, 1)

    def test_empty(self):
        np.testing.assert_array_equal(
            self.system.part.by_ids(
                []).pos, np.empty(0))

    def test_len(self):
        self.assertEqual(len(self.system.part.by_ids([])), 0)
        self.assertEqual(len(self.system.part.by_ids([0])), 1)
        self.assertEqual(len(self.system.part.by_ids([0, 1, 2])), 3)

    def test_non_existing_property(self):
        with self.assertRaises(AttributeError):
            self.all_partcls.thispropertydoesnotexist = 1.0


if __name__ == "__main__":
    ut.main()
