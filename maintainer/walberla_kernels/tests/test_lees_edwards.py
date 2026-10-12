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

import pystencils as ps
from pystencils.sympyextensions.typed_sympy import TypeCast

import lees_edwards


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


if __name__ == "__main__":
    unittest.main()
