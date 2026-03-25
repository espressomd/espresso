#
# Copyright (C) 2021-2026 The ESPResSo project
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

from pystencils_espresso import get_loop_counter_symbol

try:
    from pystencils.typing import CastFunc
except ImportError:
    from pystencils.sympyextensions import CastFunc

from pystencils import TypedSymbol
from pystencils import Assignment

from lbmpy.macroscopic_value_kernels import macroscopic_values_setter

import sympy as sp
import types


def bind_diff(self, *symbols, **kwargs):  # pylint: disable=unused-argument
    """
    Make a Symbol object behave as a constant during symbolic differentiation.
    """
    return sp.Derivative(sp.Integer(0), *symbols, **kwargs)


def bind_diff_wrt(self):  # pylint: disable=unused-argument
    """
    Make a Symbol object behave as a constant during symbolic differentiation.
    """
    return True


def velocity_offset_eqs(config, method, pdfs, shear_dir_normal, stencil,
                        with_simd, combined_kernel=True):
    """Calculates the difference between quilibrium pdf distributions
    with (rho, u) and (rho, u+v) and applies them to out-flowing
    populations in the boundary layer. Returns an AssignmentCollection
    with one Assignment per stencil direction.
    """
    dim = len(stencil[0])
    default_dtype = config.data_type.default_factory()

    # Placeholders indicating a population flows up or down.
    # Will be replaced later using the component of the stencil direction
    # along the shear_dir_normal.
    points_up = sp.Symbol('points_up')
    points_down = sp.Symbol('points_down')

    # Symbol for the coordinate index within the field,
    # used to identify boundary layers
    counters = [get_loop_counter_symbol(i) for i in range(dim)]

    lebc_top_index = TypedSymbol("lebc_top_index", counters[1].dtype)
    lebc_bot_index = TypedSymbol("lebc_bot_index", counters[1].dtype)

    # When the kernels are combined, the integration loop also reaches
    # one grid cell further to pull the populations. This changes the
    # position where to apply the LE boundary by one cell.
    boundary_offset = [0, -1]
    if combined_kernel:
        boundary_offset = [1, 0]
    if with_simd:
        # make type differentiable (CastFunc is a function taking 2 arguments,
        # but the second argument is a C++ type, i.e. a compile-time constant)
        default_dtype.diff = types.MethodType(bind_diff, default_dtype)
        default_dtype._diff_wrt = types.MethodType(
            bind_diff_wrt, default_dtype)
        # cell index conversion to floating-point is mandatory for vectorization
        # (only float32 and float64 are supported in pystencils 1.4)
        index = CastFunc(counters[1], default_dtype)
        # epsilon to avoid an off-by-one-error during floating-point comparison
        epsilon = 0.1
        index_bot = CastFunc(lebc_bot_index, default_dtype) + \
            CastFunc(boundary_offset[0] + epsilon, default_dtype)
        index_top = CastFunc(lebc_top_index, default_dtype) + \
            CastFunc(boundary_offset[1] - epsilon, default_dtype)
    else:
        index = counters[1]
        index_bot = lebc_bot_index + boundary_offset[0]
        index_top = lebc_top_index + boundary_offset[1]
    # +,-1 for upper/lower boundary layers, 0 otherwise.
    # Based on symbolic counters defined above. Only becomes
    # non-zero if the corresponding points_up/down flags
    # are engaged (which is only done for out-flowing populations)
    layer_prefactor = sp.Piecewise(
        (-1, sp.And(index <= index_bot, points_down)),
        (+1, sp.And(index >= index_top, points_up)),
        (0, True)
    )

    # Start with an equilibrium distribution for a given density and velocity
    delta_pdf_eqs = macroscopic_values_setter(
        method, sp.Symbol("dens"), [
            sp.Symbol("v_0"), sp.Symbol("v_1"), sp.Symbol("v_2")], pdfs)

    # Replace the assignments of (rho,u) by (rho, u+v) - (rho,u)
    ma = []
    for a, c in zip(delta_pdf_eqs.main_assignments, method.stencil):
        # Determine direction of the stencil component in the
        # shear_dir_normal
        if c[shear_dir_normal] == 1:
            up = True
            down = False
        elif c[shear_dir_normal] == -1:
            up = False
            down = True
        else:
            up = False
            down = False

        # Replace (rho,u) by (rho,u+v) in boundary layers
        rhs = sp.simplify(
            a.rhs -
            a.rhs.replace(
                sp.Symbol("u_0"),
                sp.Symbol("u_0") +
                layer_prefactor *
                sp.Symbol("v_s")))

        # Only engage if the population is outflowing. See layer_prefactor
        rhs = rhs.replace(points_up, up)
        rhs = rhs.replace(points_down, down)
        new_a = Assignment(a.lhs, rhs)

        ma.append(new_a)
        print(c, ma[-1])
    # Plug in modified assignments
    delta_pdf_eqs.main_assignments = ma
    return delta_pdf_eqs.main_assignments


def add_lees_edwards_to_collision(
        config, collision, pdfs, stencil, shear_dir_normal, with_simd, combined_kernel):
    # Get population shift for outflowing populations at the boundaries
    offset = velocity_offset_eqs(
        config,
        collision.method,
        pdfs,
        shear_dir_normal,
        stencil,
        with_simd,
        combined_kernel)

    ma = []
    for i, a in enumerate(collision.main_assignments):
        # Add Lees-Edwards-shift to collision main assignments
        new_a = Assignment(a.lhs, a.rhs + offset[i].rhs)
        ma.append(new_a)
    collision.main_assignments = ma
    return collision
