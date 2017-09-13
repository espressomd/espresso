from __future__ import print_function
import unittest as ut
import espressomd
from espressomd import has_features
import numpy as np


class ParticleSliceTest(ut.TestCase):

    state = [[0, 0, 0], [0, 0, 1]]
    system = espressomd.System()

    def __init__(self, *args, **kwargs):
        super(ParticleSliceTest, self).__init__(*args, **kwargs)
        self.system.part.clear()
        self.system.part.add(pos=[0, 0, 0])
        self.system.part.add(pos=[0, 0, 1])

        if has_features(["EXTERNAL_FORCES"]):
            self.system.part[1].fix = self.state[1]
            self.assertTrue(np.array_equal(
                self.system.part[0].fix, self.state[0]))
            self.assertTrue(np.array_equal(
                self.system.part[1].fix, self.state[1]))
            self.assertTrue(np.array_equal(
                self.system.part[:].fix, self.state))
        xs = self.system.part[:].pos
        for i in range(len(xs)):
            self.assertTrue(np.array_equal(xs[i], self.system.part[i].pos))

    @ut.skipIf(
        not has_features(
            ["EXTERNAL_FORCES"]),
        "Features not available, skipping test!")
    def test_1_set_different_values(self):
        self.state[0] = [1, 0, 0]
        self.state[1] = [1, 0, 0]
        self.system.part[:].fix = self.state
        self.assertTrue(np.array_equal(self.system.part[:].fix, self.state))

    @ut.skipIf(
        not has_features(
            ["EXTERNAL_FORCES"]),
        "Features not available, skipping test!")
    def test_2_set_same_value(self):
        self.state[0] = [0, 1, 0]
        self.state[1] = [0, 1, 0]
        self.system.part[:].fix = self.state[1]
        self.assertTrue(np.array_equal(self.system.part[:].fix, self.state))

    @ut.skipIf(
        not has_features(
            ["EXTERNAL_FORCES"]),
        "Features not available, skipping test!")
    def test_3_set_one_value(self):
        self.state[1] = [0, 0, 1]
        self.system.part[1:].fix = self.state[1]
        self.assertTrue(np.array_equal(self.system.part[:].fix, self.state))

    @ut.skipIf(
        not has_features(
            ["EXTERNAL_FORCES"]),
        "Features not available, skipping test!")
    def test_4_str(self):
        self.assertEqual(
            repr(self.system.part[0].fix), repr(np.array([0, 1, 0])))
        self.assertEqual(repr(self.system.part[:].fix), repr(
            np.array([[0, 1, 0], [0, 0, 1]])))

    def test_pos_str(self):
        self.assertEqual(repr(self.system.part[0].pos), repr(
            np.array([0.0, 0.0, 0.0])))
        if (len(self.system.part)) > 1:
            self.assertEqual(repr(self.system.part[:].pos), repr(
                np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]])))

    @ut.skipIf(
        not has_features(
            ["ELECTROSTATICS"]),
        "Features not available, skipping test!")
    def test_scalar(self):
        self.system.part[:1].q = 1.3
        self.assertEqual(self.system.part[0].q, 1.3)
        self.system.part[:].q = 2.0
        self.assertEqual(self.system.part[0].q, 2)
        self.assertEqual(self.system.part[1].q, 2)
        self.system.part[:].q = 3
        self.assertEqual(self.system.part[0].q, 3)
        self.assertEqual(self.system.part[1].q, 3)
        self.system.part[:].q = [-1, 1.0]
        self.assertEqual(self.system.part[0].q, -1)
        self.assertEqual(self.system.part[1].q, 1)
        qs = self.system.part[:].q
        self.assertEqual(qs[0], -1)
        self.assertEqual(qs[1], 1)

    def test_empty(self):
        self.assertTrue(np.array_equal(self.system.part[0:0].pos, np.empty(0)))


if __name__ == "__main__":
    ut.main()
