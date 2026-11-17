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

import unittest

import sympy as sp
import pystencils as ps
import lbmpy
import lbmpy.creationfunctions
import lbmpy.forcemodels
from pystencils.sympyextensions.typed_sympy import TypeCast

import lees_edwards
import pystencils_espresso
import relaxation_rates


class BoundaryLayerMasksTest(unittest.TestCase):

    def evaluate(self, mask, index, bot, top):
        values = {"ctr_1": index, "lebc_bot_index": bot,
                  "lebc_top_index": top}
        expr = mask.replace(lambda e: isinstance(e, TypeCast),
                            lambda e: e.args[0])
        return int(expr.subs({s: values[s.name] for s in expr.free_symbols}))

    def test_masks_select_the_outermost_layers(self):
        is_top, is_bot = lees_edwards.boundary_layer_masks(
            1, ps.create_type("float64"))
        bot, top = 0, 10
        for index in range(-1, top + 1):
            self.assertEqual(self.evaluate(is_top, index, bot, top),
                             int(index >= top - 1), index)
            self.assertEqual(self.evaluate(is_bot, index, bot, top),
                             int(index <= bot), index)

    def test_masks_on_a_block_inside_the_box(self):
        is_top, is_bot = lees_edwards.boundary_layer_masks(
            1, ps.create_type("float64"))
        # block covering the cell layers 8..15 of a box of thickness 16
        bot, top = 0 - 8, 16 - 8
        rows = range(8)
        self.assertEqual([self.evaluate(is_top, i, bot, top) for i in rows],
                         [0] * 7 + [1])
        self.assertEqual([self.evaluate(is_bot, i, bot, top) for i in rows],
                         [0] * 8)


class AddLeesEdwardsToCollisionTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        stencil = lbmpy.LBStencil(lbmpy.Stencil.D3Q19)
        dtype = ps.create_type("float64")
        fields = pystencils_espresso.generate_fields(stencil, dtype)
        method = lbmpy.creationfunctions.create_mrt_orthogonal(
            stencil=stencil, relaxation_rates=relaxation_rates.rr_getter,
            force_model=lbmpy.forcemodels.Schiller(
                fields["force"].center_vector),
            weighted=True, compressible=True, delta_equilibrium=True,
            zero_centered=True)
        config = lbmpy.LBMConfig(
            stencil=stencil, streaming_pattern="pull", weighted=True,
            compressible=True, delta_equilibrium=True, zero_centered=True)
        opt = lbmpy.LBMOptimisation(
            symbolic_field=fields["pdfs"],
            symbolic_temporary_field=fields["pdfs_tmp"], field_layout="fzyx")
        collision = lbmpy.creationfunctions.create_lb_collision_rule(
            method, lbm_config=config, lbm_optimisation=opt)
        cls.original = list(collision.main_assignments)
        cls.collision = lees_edwards.add_lees_edwards_to_collision(
            collision, fields["pdfs"], 1, dtype)
        cls.stencil = method.stencil
        cls.deltas = [new.rhs - old.rhs for new, old in
                      zip(collision.main_assignments, cls.original)]
        cls.le_top, cls.le_bot = sp.symbols("le_top le_bot")

    def shifted(self, deltas):
        return {i for i, d in enumerate(deltas) if sp.simplify(d) != 0}

    def test_shifted_directions(self):
        # the equilibrium of the populations with c_y = 0, or with only c_y
        # non-zero and c_x = 0, doesn't depend on u_0 at second order
        expected = {i for i, c in enumerate(self.stencil)
                    if c[1] != 0 and (c[0] != 0 or c[2] == 0)}
        self.assertEqual(self.shifted(self.deltas), expected)

    def test_top_and_bottom_are_independent(self):
        all_shifted = self.shifted(self.deltas)
        up = {i for i in all_shifted if self.stencil[i][1] == 1}
        down = {i for i in all_shifted if self.stencil[i][1] == -1}
        self.assertTrue(up and down)
        no_top = [d.subs(self.le_top, 0) for d in self.deltas]
        no_bot = [d.subs(self.le_bot, 0) for d in self.deltas]
        self.assertEqual(self.shifted(no_top), down)
        self.assertEqual(self.shifted(no_bot), up)

    def test_shift_uses_collision_symbols(self):
        symbols = set()
        for d in self.deltas:
            symbols |= d.free_symbols
        symbols -= {sp.Symbol("v_s"), self.le_top, self.le_bot}
        known = {a.lhs for a in self.collision.subexpressions}
        self.assertTrue(symbols)
        self.assertTrue(symbols <= known, symbols - known)

    def test_mask_subexpressions_come_first(self):
        self.assertEqual(self.collision.subexpressions[0].lhs, self.le_top)
        self.assertEqual(self.collision.subexpressions[1].lhs, self.le_bot)


if __name__ == "__main__":
    unittest.main()
