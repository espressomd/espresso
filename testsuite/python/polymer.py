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
from espressomd.shapes import Wall


class PolymerPositions(ut.TestCase):
    box_l = 15
    seed = 23

    system = espressomd.System(box_l=[box_l, box_l, box_l])
    np.random.seed(1234)
    system.set_random_state_PRNG()

    def assertShape(self, positions, n_poly, n_mono):
        """
        Assert that positions array has expected shape and
        expected number of elements.

        """
        shape = np.shape(positions)
        self.assertEqual(shape, (n_poly, n_mono, 3))

    def assertBondLength(self, positions, bond_length):
        """
        Assert that distance between monomers is bond_length.

        """
        for p in positions:
            distances = np.linalg.norm(p[1:] - p[:-1], axis=1)
            max_deviation = np.max(abs(distances - bond_length))
            self.assertLess(max_deviation, 1e-10)

    def assertBondAngle(self, positions, bond_angle, bond_length):
        """
        Assert that angle between adjacent beads in positions is bond_angle.

        """
        for p in positions:
            distance_vectors = (p[1:] - p[:-1]) / bond_length
            cos_angles = np.einsum(
                'ij, ij -> i',
                distance_vectors[1:],
                -distance_vectors[:-1])
            max_deviation = max(abs(cos_angles - np.cos(bond_angle)))
            self.assertLess(max_deviation, 1e-10)

    def assertMinDistGreaterEqual(self, positions, r_min):
        """
        Assert that all points in positions are at least r_min away from each other.

        """
        particle_positions = positions.flatten().reshape((-1, 3))
        # use folded coordinates
        particle_positions %= self.box_l
        for pos in particle_positions:
            distances = np.linalg.norm(particle_positions - pos, axis=1)
            # exclude zero distance, i.e.,  distance to pos itself
            distances = distances[np.nonzero(distances)]
            self.assertGreaterEqual(min(distances), r_min)

    def test_bond_lengths(self):
        """
        Check that distance between neighboring monomers is indeed bond_length.

        """
        bond_lengths = [0.735, 1.459]
        num_poly = 10
        num_mono = 25
        for bond_length in bond_lengths:
            positions = polymer.positions(
                n_polymers=num_poly, beads_per_chain=num_mono,
                    bond_length=bond_length, seed=self.seed)

            self.assertShape(positions, num_poly, num_mono)
            self.assertBondLength(positions, bond_length)

    def test_bond_angles(self):
        """
        Check that bond_angle is obeyed.

        """
        bond_angles = [0.436 * np.pi, np.pi / 3., np.pi / 5.]
        num_poly = 10
        num_mono = 25
        bond_length = 1.34
        for bond_angle in bond_angles:
            positions = polymer.positions(
                n_polymers=num_poly, beads_per_chain=num_mono,
                    bond_angle=bond_angle,
                    bond_length=bond_length, seed=self.seed)

            self.assertShape(positions, num_poly, num_mono)
            self.assertBondLength(positions, bond_length)
            self.assertBondAngle(positions, bond_angle, bond_length)

    def test_start_positions(self):
        """
        Check that setting start positions behaves correctly.

        """
        num_poly = 90
        num_mono = 25
        bond_length = 0.83
        start_positions = np.random.random((num_poly, 3)) * self.box_l

        # make sure that incorrect size leads to error
        with self.assertRaises(ValueError):
            positions = polymer.positions(
                n_polymers=num_poly + 1,
                    beads_per_chain=num_mono,
                    start_positions=start_positions,
                    bond_length=bond_length, seed=self.seed)

        # check that start positions are actually used
        positions = polymer.positions(
            n_polymers=num_poly,
                beads_per_chain=num_mono,
                start_positions=start_positions,
                bond_length=bond_length, seed=self.seed)

        self.assertListEqual(
            start_positions.tolist(),
            positions[:,
                      0].tolist())

    def test_min_dist(self):
        """
        Check that min_dist is respected.

        """
        num_poly = 5
        num_mono = 150
        bond_length = 0.945

        positions = polymer.positions(
            n_polymers=num_poly,
                beads_per_chain=num_mono,
                bond_length=bond_length,
                min_distance=bond_length, seed=self.seed)

        self.assertBondLength(positions, bond_length)
        self.assertMinDistGreaterEqual(positions, bond_length - 1e-10)

    def test_respect_constraints_wall(self):
        """
        Check that constraints are respected.

        """
        num_poly = 20
        num_mono = 5
        bond_length = 1.19

        w = espressomd.shapes.Wall(normal=[0., 0., 1.], dist=0.5 * self.box_l)
        wall_constraint = espressomd.constraints.ShapeBasedConstraint(shape=w)
        c = self.system.constraints.add(wall_constraint)

        positions = polymer.positions(
            n_polymers=num_poly,
                beads_per_chain=num_mono,
                bond_length=bond_length,
                respect_constraints=True, seed=self.seed)

        positions %= self.box_l

        z_components = positions[:, :, 2][0]
        for z in z_components:
            print(z)
            self.assertGreaterEqual(z, 0.5 * self.box_l)

        # assert that illegal start position raises error
        assertRaisesRegex = self.assertRaisesRegexp if sys.version_info < (
            3, 2) else self.assertRaisesRegex
        with assertRaisesRegex(Exception, 'Invalid start positions.'):
            illegal_start = np.array([[1., 1., 0.2 * self.box_l]])
            positions = polymer.positions(
                n_polymers=1,
                    beads_per_chain=10,
                    start_positions=illegal_start,
                    bond_length=bond_length,
                    respect_constraints=True, seed=self.seed)
        self.system.constraints.remove(wall_constraint)


if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
