#
# Copyright (C) 2026 The ESPResSo project
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

"""Tests of the layout of the flux field of the electrokinetics."""

import itertools
import unittest

import pystencils as ps
import sympy as sp

import flux_layout


class FluxLayoutTest(unittest.TestCase):

    def test_sizes(self):
        self.assertEqual(flux_layout.DIM, 3)
        self.assertEqual(flux_layout.FLUX_COUNT, 13)
        self.assertEqual(len(flux_layout.STORED_DIRECTIONS), 13)
        self.assertEqual(len(flux_layout.NEIGHBORS), 26)

    def test_neighbors_are_the_directions_to_the_26_neighbors(self):
        expected = {d for d in itertools.product((-1, 0, 1), repeat=3)
                    if any(d)}
        self.assertEqual(set(flux_layout.NEIGHBORS), expected)

    def test_a_cell_stores_one_direction_of_each_pair(self):
        stored = set(flux_layout.STORED_DIRECTIONS)
        self.assertEqual(len(stored), 13)
        inverse = {tuple(-c for c in d) for d in stored}
        self.assertEqual(stored & inverse, set())
        self.assertEqual(stored | inverse, set(flux_layout.NEIGHBORS))

    def test_layout_is_the_layout_of_pystencils(self):
        # the components of the flux field of the bridge (13 components of a
        # GhostLayerField) have this order; the layout of pystencils 1.4
        # (negative directions, another order) is not compatible
        self.assertEqual(flux_layout.STORED_DIRECTIONS, (
            (0, 1, 0), (1, 0, 0), (0, 0, 1), (1, 1, 0), (1, -1, 0),
            (0, 1, 1), (1, 0, 1), (0, 1, -1), (1, 0, -1), (1, 1, 1),
            (-1, 1, 1), (1, -1, 1), (1, 1, -1)))

    def test_flux_field(self):
        field = flux_layout.flux_field(ps.create_type("float32"), "j")
        self.assertEqual(field.name, "j")
        self.assertEqual(field.spatial_dimensions, 3)
        self.assertEqual(field.index_shape, (13,))
        self.assertEqual(field.layout, (2, 1, 0))
        self.assertEqual(field.dtype, ps.create_type("float32"))

    def test_face_of_a_stored_direction(self):
        for k, direction in enumerate(flux_layout.STORED_DIRECTIONS):
            self.assertEqual(flux_layout.face_location(direction),
                             ((0, 0, 0), k, 1))

    def test_face_of_the_other_directions_is_in_the_neighbor(self):
        for direction in flux_layout.NEIGHBORS:
            if direction in flux_layout.STORED_DIRECTIONS:
                continue
            offset, component, sign = flux_layout.face_location(direction)
            self.assertEqual(offset, direction)
            self.assertEqual(sign, -1)
            self.assertEqual(
                flux_layout.STORED_DIRECTIONS[component],
                tuple(-c for c in direction))

    def test_both_sides_of_a_face_are_the_same_element(self):
        # the face between c and c + d seen from c (direction d) and from
        # c + d (direction -d): the same component of the same cell, with
        # opposite signs (the flux along d is the opposite of the one along
        # -d)
        for direction in flux_layout.NEIGHBORS:
            offset, component, sign = flux_layout.face_location(direction)
            owner = offset
            inverse = tuple(-c for c in direction)
            offset2, component2, sign2 = flux_layout.face_location(inverse)
            owner2 = tuple(d + o for d, o in zip(direction, offset2))
            self.assertEqual(owner, owner2, direction)
            self.assertEqual(component, component2, direction)
            self.assertEqual(sign, -sign2, direction)

    def test_every_face_of_a_cell_is_a_different_element(self):
        elements = set()
        for direction in flux_layout.NEIGHBORS:
            offset, component, _ = flux_layout.face_location(direction)
            elements.add((offset, component))
        self.assertEqual(len(elements), 26)

    def test_no_face_of_the_cell_itself_or_far_away(self):
        for direction in ((0, 0, 0), (2, 0, 0), (1, 1, 2)):
            with self.assertRaises(ValueError):
                flux_layout.face_location(direction)

    def test_lower_replaces_the_staggered_accesses(self):
        dtype = ps.create_type("float64")
        plain = flux_layout.flux_field(dtype, "j")
        staggered = flux_layout.staggered_flux_field(dtype, "j")
        expression = 2 * staggered.face[0, 0, 0, (-1, 0, 0)] \
            + staggered.face[0, 0, 0, (0, 1, 0)] + sp.Symbol("a")
        lowered = flux_layout.lower(expression, plain)
        accesses = lowered.atoms(ps.Field.Access)
        self.assertEqual(
            {(a.offsets, a.index[0]) for a in accesses},
            {((-1, 0, 0), 1), ((0, 0, 0), 0)})
        self.assertEqual({a.field for a in accesses}, {plain})
        # the opposite direction carries the opposite sign
        by_offset = {a.offsets: a for a in accesses}
        self.assertEqual(lowered.coeff(by_offset[(-1, 0, 0)]), -2)
        self.assertEqual(lowered.coeff(by_offset[(0, 0, 0)]), 1)

    def test_lowered_faces_agree_with_face_location(self):
        dtype = ps.create_type("float64")
        plain = flux_layout.flux_field(dtype, "j")
        staggered = flux_layout.staggered_flux_field(dtype, "j")
        for direction in flux_layout.NEIGHBORS:
            lowered = flux_layout.lower(
                staggered.face[0, 0, 0, direction], plain)
            (access,) = lowered.atoms(ps.Field.Access)
            offset, component, sign = flux_layout.face_location(direction)
            self.assertEqual(access.offsets, offset)
            self.assertEqual(access.index[0], component)
            self.assertEqual(lowered.coeff(access), sign)


if __name__ == "__main__":
    unittest.main()
