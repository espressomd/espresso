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
Reaction kernels of the electrokinetics: the sweeps
``ReactionKernelBulk<N>`` and ``ReactionKernelIndexed<N>`` of
``walberla::pystencils`` for ``N`` reactants (the species of a reaction).

A reaction changes the densities of its reactants in a cell by

    rho_i += stoech_i * rate_coefficient * prod_j rho_j ** order_j

where the products use the densities before the reaction. The kernel reads
all densities into local variables first, so that all reactants see the
same densities of the cell.

* The bulk sweep runs on all cells of the interior of a block.
* The indexed sweep runs on the cells of a
  ``walberla::v8::sweep::SparseIndexList`` that its constructor takes, which
  holds the cells with a reaction (see `sfg_driver.generate_sweep`): the
  density fields are accessed with the absolute coordinates of these cells.
  The list is in unified memory on CUDA.

The properties of the sweep are the density fields ``rho_<i>`` (BlockDataIDs)
and the scalars ``order<i>``, ``rate_coefficient`` and ``stoech<i>``, in the
data type of the kernel. The scalars are named ``order0``, ``stoech0``, ...
and not ``order_0``, ``stoech_0``, ...: sweepgen combines properties
``<name>_0`` to ``<name>_2`` to a ``Vector3``, and ``<name>_0`` and
``<name>_1`` to a ``Vector2``, which would give the classes of 2 and 3
reactants a different constructor than the classes of 1, 4 and 5 reactants.
The constructor takes the index list (indexed sweep), the BlockDataIDs of the
densities, then the scalars, each group in alphabetical order, like the
constructors of the old ``pystencils_walberla`` classes did.
"""

from pystencils import Assignment, AssignmentCollection, Field, TypedSymbol
from pystencils import FieldType
import sympy as sp

import sfg_driver

MAX_NUM_REACTANTS = 5


def reaction_assignments(dtype, num_reactants, indexed=False):
    """
    The assignments of a reaction of ``num_reactants`` reactants, of the
    data type ``dtype``. With ``indexed``, the cell is given by the index
    list of the sparse sweep (`sfg_driver.index_list_field`) and the density
    fields have the type ``FieldType.CUSTOM``, which are accessed with
    absolute coordinates.
    """
    if not 1 <= num_reactants <= MAX_NUM_REACTANTS:
        raise ValueError(
            f"reactions of {num_reactants} reactants are not implemented, "
            f"1 to {MAX_NUM_REACTANTS} are")
    kind = FieldType.CUSTOM if indexed else FieldType.GENERIC
    rhos = [Field.create_generic(f"rho_{i}", 3, dtype, index_dimensions=0,
                                 layout="fzyx", field_type=kind)
            for i in range(num_reactants)]
    if indexed:
        cell = tuple(sfg_driver.index_list_field()[0](c) for c in "xyz")
        cells = [rho.absolute_access(cell, ()) for rho in rhos]
    else:
        cells = [rho.center for rho in rhos]
    orders = [TypedSymbol(f"order{i}", dtype) for i in range(num_reactants)]
    stoech = [TypedSymbol(f"stoech{i}", dtype) for i in range(num_reactants)]
    rate_coefficient = TypedSymbol("rate_coefficient", dtype)
    local_rhos = sp.symbols(f"local_rho_:{num_reactants}")
    rate_factor = sp.Symbol("rate_factor")
    rate = rate_coefficient
    for local_rho, order in zip(local_rhos, orders):
        rate *= sp.Pow(local_rho, order)
    subexpressions = [Assignment(local_rho, cell)
                      for local_rho, cell in zip(local_rhos, cells)]
    subexpressions.append(Assignment(rate_factor, rate))
    main_assignments = [Assignment(cell, local_rho + rate_factor * coef)
                        for cell, local_rho, coef in
                        zip(cells, local_rhos, stoech)]
    return AssignmentCollection(main_assignments, subexpressions)


def generate_bulk(sfg, ctx, name, num_reactants):
    """Generate the bulk reaction sweep ``name`` (dense, on the interior)."""
    sfg_driver.generate_sweep(
        sfg, name, reaction_assignments(ctx.data_type, num_reactants))


def generate_indexed(sfg, ctx, name, num_reactants):
    """Generate the indexed reaction sweep ``name`` (sparse, index list)."""
    sfg_driver.generate_sweep(
        sfg, name,
        reaction_assignments(ctx.data_type, num_reactants, indexed=True),
        sparse=True)
