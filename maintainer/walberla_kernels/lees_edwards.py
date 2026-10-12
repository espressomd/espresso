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
Lees-Edwards boundary conditions for the LB stream-collide kernel.

In the outermost cell layers of the shear-normal direction, the populations
that flow out of the simulation box are shifted by the difference between
the equilibria with velocity ``u`` and ``u + v_s``, where ``v_s`` is the
shear velocity.
"""

import sympy as sp
from pystencils import Assignment, DynamicType, TypedSymbol
from pystencils.defaults import DEFAULTS
from pystencils.sympyextensions import tcast
from lbmpy.macroscopic_value_kernels import macroscopic_values_setter


def boundary_layer_masks(shear_dir_normal, dtype):
    """
    Masks that are 1 in the top and bottom cell layer, respectively, and 0
    elsewhere. ``lebc_bot_index`` and ``lebc_top_index`` are the block-local
    indices of the first cell of the simulation box and of the first cell
    after it along the shear-normal direction.

    The masks are written as ``max(0, min(1, ...))`` of the integer cell
    index instead of a ``Piecewise``: they only depend on the cell layer, so
    pystencils hoists them out of the vectorized innermost loop, while a
    ``Piecewise`` can't be vectorized by pystencils 2.1.
    """
    index = tcast(DEFAULTS.spatial_counters[shear_dir_normal], dtype)
    top = tcast(TypedSymbol("lebc_top_index", DynamicType.INDEX_TYPE), dtype)
    bot = tcast(TypedSymbol("lebc_bot_index", DynamicType.INDEX_TYPE), dtype)
    is_top = sp.Max(0, sp.Min(1, index - top + 2))  # index >= top - 1
    is_bot = sp.Max(0, sp.Min(1, bot - index + 1))  # index <= bot
    return is_top, is_bot


def add_lees_edwards_to_collision(collision, pdfs, shear_dir_normal, dtype):
    """
    Add the Lees-Edwards velocity shift of the outflowing populations to the
    main assignments of ``collision``. The shear velocity is the kernel
    parameter ``v_s``, the shear direction is x.
    """
    method = collision.method
    is_top, is_bot = boundary_layer_masks(shear_dir_normal, dtype)
    le_top, le_bot = sp.symbols("le_top le_bot")
    velocity = [sp.Symbol(f"v_{i}") for i in range(method.dim)]
    equilibrium = macroscopic_values_setter(
        method, sp.Symbol("dens"), velocity, pdfs)
    v_s = sp.Symbol("v_s")
    u_0 = sp.Symbol("u_0")
    main_assignments = []
    for a, eq, c in zip(collision.main_assignments,
                        equilibrium.main_assignments, method.stencil):
        if c[shear_dir_normal] == 1:
            shift = le_top * v_s
        elif c[shear_dir_normal] == -1:
            shift = -le_bot * v_s
        else:
            shift = 0
        delta = sp.simplify(eq.rhs - eq.rhs.replace(u_0, u_0 + shift))
        main_assignments.append(Assignment(a.lhs, a.rhs + delta))
    collision.main_assignments = main_assignments
    collision.subexpressions = [Assignment(le_top, is_top),
                                Assignment(le_bot, is_bot),
                                *collision.subexpressions]
    return collision
