#
# Copyright (C) 2022-2026 The ESPResSo project
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

"""
Layout of the flux field of the electrokinetics.

The flux field has 13 components for each cell (the memory layout of a
``walberla::field::GhostLayerField<T, 13>`` in ``fzyx`` order). The component
``k`` of the cell ``c`` is the flux through the face between ``c`` and
``c + e_k`` along ``e_k``, where ``e_k`` is the k-th direction of
`STORED_DIRECTIONS`: half of the 26 neighbors of a cell (the *positive*
edges, one of each pair of opposite directions). The flux through the other
faces of a cell is stored in its neighbors, with the opposite sign. This is
the layout of ``Staggering.Corner3D`` of pystencils 2.1, which differs from
the layout of the staggered fields of pystencils 1.4 (the *negative* edges,
in another order): a flux field is not compatible with the one of earlier
versions.

pystencils 2.1 derives the finite-volume expressions from a symbolic
`pystencils.StaggeredField`. sweepgen can't handle staggered fields, so the
accesses to the staggered field are lowered to accesses to a plain `Field`
with 13 components (`lower`), which is the type that the bridge allocates.

This module is the single place that knows the layout: the sweeps
(`ek_sweeps`), the fixed flux boundary (`boundaries`) and the accessors
(`ek_accessors`) use it, and all of them must change together if the layout
changes.
"""

from pystencils import Field
from pystencils.grids import StaggeredField, StaggeredFieldAccess, Staggering

DIM = 3
STAGGERING = Staggering.Corner3D
# the direction of each component of the flux field
STORED_DIRECTIONS = tuple(STAGGERING.staggered_entries)
FLUX_COUNT = len(STORED_DIRECTIONS)
# the directions to the 26 neighbors of a cell
NEIGHBORS = tuple(STAGGERING.stencil_entries[1:])


def flux_field(dtype, name="flux"):
    """The plain field with the 13 components of the flux of a cell."""
    return Field.create_generic(
        name, DIM, dtype, index_dimensions=1, index_shape=(FLUX_COUNT,),
        layout="fzyx")


def staggered_flux_field(dtype, name="flux"):
    """
    The symbolic staggered flux field, which the finite-volume derivation of
    pystencils needs. Its accesses must be lowered with `lower` to the
    accesses of `flux_field` of the same name.
    """
    return StaggeredField(name, STAGGERING, (), dtype=dtype)


def lower(expression, plain_field):
    """
    Replace the accesses to a staggered flux field in ``expression`` by the
    accesses to ``plain_field`` (see `flux_field`) of the same cell and
    component.
    """
    def is_access(e):
        return isinstance(e, StaggeredFieldAccess)

    def lower_access(e):
        assert e.vector_index is None
        offsets = tuple(int(o) for o in e.offsets)
        return plain_field[offsets](int(e.index))

    return expression.replace(is_access, lower_access)


def face_location(direction):
    """
    Where the flux through the face between a cell ``c`` and its neighbor
    ``c + direction`` along ``direction`` is stored: returns the tuple
    ``(offset, component, sign)`` of the cell ``c + offset`` that holds it,
    the component of the flux field, and the sign, such that the flux is
    ``sign * flux(c + offset)[component]``.

    :raises ValueError: if ``direction`` isn't a direction to a neighbor
    """
    direction = tuple(int(d) for d in direction)
    if direction not in NEIGHBORS:
        raise ValueError(f"{direction} is not a direction to a neighbor")
    if direction in STORED_DIRECTIONS:
        return (0, 0, 0), STORED_DIRECTIONS.index(direction), 1
    opposite = tuple(-d for d in direction)
    return direction, STORED_DIRECTIONS.index(opposite), -1
