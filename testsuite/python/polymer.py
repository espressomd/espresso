# Copyright (C) 2010-2018 The ESPResSo project
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
from __future__ import print_function
import sys
import unittest as ut
import numpy as np
import espressomd
from espressomd.interactions import FeneBond
from espressomd import polymer

class PolymerPositions(ut.TestCase):
    box_l = 15

    system = espressomd.System(box_l=[box_l, box_l, box_l])
    np.random.seed(1234)
    system.set_random_state_PRNG()

    def assertShape(self, positions, n_poly, n_mono):
        """Assert that positions array has expected shape and
        expected number of elements.

        """
        shape = np.shape(positions)
        self.assertEqual(shape, (n_poly, n_mono, 3))

    def assertBondLength(self, positions, bond_length):
        for p in positions:
            distances = np.linalg.norm(p[1:] - p[:-1], axis=1)
            max_deviation = np.max(abs(distances - bond_length))
            self.assertLess(max_deviation, 1e-10)

    def assertBondAngle(self, positions, bond_angle, bond_length):
        for p in positions:
            distance_vectors = (p[1:] - p[:-1]) / bond_length
            cos_angles = np.einsum('ij, ij -> i', distance_vectors[1:], distance_vectors[:-1])
            max_deviation = max(abs(cos_angles - np.cos(bond_angle)))
            self.assertLess(max_deviation, 1e-2)

    def test_bond_lengths(self):
        """Check that distance between neighboring monomers is indeed bond_length.

        """
        bond_lengths = [0.735, 1.459]
        num_poly = 10
        num_mono = 25
        for bond_length in bond_lengths:
            positions = polymer.polymer_positions(
                    polymers=num_poly, monomers=num_mono,
                    bond_length=bond_length,
                    max_tries=3000)

            self.assertShape(positions, num_poly, num_mono)
            self.assertBondLength(positions, bond_length)

    def test_bond_angles(self):
        """Check that bond_angle is obeyed.

        """
        bond_angles = [0.436 * np.pi, np.pi/3., np.pi/5.]
        num_poly = 10
        num_mono = 25
        bond_length = 1.34
        for bond_angle in bond_angles:
            positions = polymer.polymer_positions(
                    polymers=num_poly, monomers=num_mono,
                    bond_angle=bond_angle,
                    bond_length=bond_length,
                    max_tries=15000)

            self.assertShape(positions, num_poly, num_mono)
            self.assertBondLength(positions, bond_length)
            self.assertBondAngle(positions, bond_angle, bond_length)

    def test_start_positions(self):
        """Check that setting start positions behaves correctly.

        """
        num_poly = 10
        num_mono = 25
        bond_length = 0.83
        start_positions = np.random.random((num_poly, 3)) * self.box_l

        # make sure that incorrect size leads to error
        with self.assertRaises(ValueError):
            positions = polymer.polymer_positions(
                    polymers=num_poly+1,
                    monomers=num_mono,
                    start_positions=start_positions,
                    bond_length=bond_length)

        # check that start positions are actually used
        positions = polymer.polymer_positions(
                polymers=num_poly,
                monomers=num_mono,
                start_positions=start_positions,
                bond_length=bond_length)

        self.assertListEqual(start_positions.tolist(), positions[:,0].tolist())


if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
