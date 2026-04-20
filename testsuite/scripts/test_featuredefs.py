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

import io
import sys
import graphlib
import importlib
import unittest as ut

sys.path.insert(0, "@CMAKE_SOURCE_DIR@/src/config")
module = importlib.import_module("featuredefs")


class Test(ut.TestCase):

    def test_parser(self):
        ruleset = "A\nB implies C,D,\nE requires F and not G\nG equals H or J and not K\nL external\nM external"
        obj = module.defs(io.StringIO(ruleset))
        self.assertEqual(obj.allfeatures, set("ABEGLM"))
        self.assertEqual(obj.implications, [("B", "C,D")])
        self.assertEqual(
            obj.requirements, [
                ("E", "F and not G", "defined(ESPRESSO_F) && ! defined(ESPRESSO_G)")])
        self.assertEqual(obj.derived, set("G"))
        self.assertEqual(
            obj.derivations, [
                ("G", "H or J and not K",
                 "defined(ESPRESSO_H) || defined(ESPRESSO_J) && ! defined(ESPRESSO_K)")])
        self.assertEqual(obj.externals, set("LM"))

    def test_implications_topological_ordering(self):
        # ruleset: A -> B
        # round 1
        # B must be resolved by A
        ruleset = "A implies B"
        round1 = {"B": {"A"}}
        obj = module.defs(io.StringIO(ruleset))
        self.assertEqual(obj.process_implications(), [round1])
        # ruleset: A -> B, C -> D
        # round 1
        # B must be resolved by A
        # D must be resolved by C
        ruleset = "A implies B\nC implies D"
        round1 = {"B": {"A"}, "D": {"C"}}
        obj = module.defs(io.StringIO(ruleset))
        self.assertEqual(obj.process_implications(), [round1])
        # ruleset: A -> C, B -> C
        # round 1
        # C must be resolved by A
        # C must be resolved by B
        ruleset = "A implies C\nB implies C"
        round1 = {"C": {"A", "B"}}
        obj = module.defs(io.StringIO(ruleset))
        self.assertEqual(obj.process_implications(), [round1])
        # ruleset: A -> B, B -> C
        # round 1
        # B must be resolved by A
        # round 2
        # C must be resolved by B
        ruleset = "A implies B\nB implies C"
        round1 = {"B": {"A"}}
        round2 = {"C": {"B"}}
        obj = module.defs(io.StringIO(ruleset))
        self.assertEqual(obj.process_implications(), [round1, round2])
        # ruleset: A -> B, B -> {C, D}, E -> C
        # round 1
        # B must be resolved by A
        # C must be resolved by E
        # round 2
        # D must be resolved by B
        # C must be resolved by B
        ruleset = "A implies B\nB implies C, D\nE implies C"
        round1 = {"B": {"A"}, "C": {"E"}}
        round2 = {"C": {"B"}, "D": {"B"}}
        obj = module.defs(io.StringIO(ruleset))
        self.assertEqual(obj.process_implications(), [round1, round2])
        # ruleset: A -> B, B -> {C, D}, E -> C, F -> C
        # round 1
        # B must be resolved by A
        # C must be resolved by E
        # C must be resolved by F
        # round 2
        # D must be resolved by B
        # C must be resolved by B
        ruleset = "A implies B\nB implies C, D\nE implies C\nF implies C"
        round1 = {"B": {"A"}, "C": {"E", "F"}}
        round2 = {"C": {"B"}, "D": {"B"}}
        obj = module.defs(io.StringIO(ruleset))
        self.assertEqual(obj.process_implications(), [round1, round2])
        # ruleset: A -> {B, G}, B -> {C, D}, E -> C, F -> C, G -> C
        # round 1
        # B must be resolved by A
        # G must be resolved by A
        # C must be resolved by E
        # C must be resolved by F
        # round 2
        # D must be resolved by B
        # C must be resolved by B
        # C must be resolved by G
        ruleset = "A implies B, G\nB implies C, D\nE implies C\nF implies C\nG implies C"
        round1 = {"B": {"A"}, "C": {"E", "F"}, "G": {"A"}}
        round2 = {"C": {"B", "G"}, "D": {"B"}}
        obj = module.defs(io.StringIO(ruleset))
        self.assertEqual(obj.process_implications(), [round1, round2])
        # ruleset: A -> {B, C}, B -> {C, D}, E -> C, F -> C
        # round 1
        # B must be resolved by A
        # C must be resolved by A
        # C must be resolved by E
        # C must be resolved by F
        # round 2
        # D must be resolved by B
        # C must be resolved by B
        ruleset = "A implies B, C\nB implies C, D\nE implies C\nF implies C"
        round1 = {"B": {"A"}, "C": {"E", "F", "A"}}
        round2 = {"C": {"B"}, "D": {"B"}}
        obj = module.defs(io.StringIO(ruleset))
        self.assertEqual(obj.process_implications(), [round1, round2])
        # ruleset: A -> {B, C}, B -> D, C -> D, D -> E
        # round 1
        # B must be resolved by A
        # C must be resolved by A
        # round 2
        # D must be resolved by B
        # D must be resolved by C
        # round 3
        # E must be resolved by D
        ruleset = "A implies B, C\nB implies D\nC implies D\nD implies E"
        round1 = {"B": {"A"}, "C": {"A"}}
        round2 = {"D": {"B", "C"}}
        round3 = {"E": {"D"}}
        obj = module.defs(io.StringIO(ruleset))
        self.assertEqual(obj.process_implications(), [round1, round2, round3])

    def test_exceptions(self):
        # detect syntax errors
        with self.assertRaisesRegex(module.DefinitionError, "<io.StringIO>:  0: <feature> equals <expr> in the following line:\nA equals"):
            module.defs(io.StringIO("A equals "))
        with self.assertRaisesRegex(module.DefinitionError, "<io.StringIO>:  1: Derived feature is already defined above: in the following line:\nA equals C"):
            module.defs(io.StringIO("A equals B\nA equals C"))
        with self.assertRaisesRegex(module.DefinitionError, "<io.StringIO>:  1: Derived feature is already defined as external above: in the following line:\nA equals B"):
            module.defs(io.StringIO("A external\nA equals B"))
        with self.assertRaisesRegex(module.DefinitionError, "<io.StringIO>:  0: <feature> external in the following line:\nA external B"):
            module.defs(io.StringIO("A external B"))
        with self.assertRaisesRegex(module.DefinitionError, "<io.StringIO>:  1: External feature is already defined as derived above: in the following line:\nA external"):
            module.defs(io.StringIO("A equals B\nA external"))
        with self.assertRaisesRegex(module.DefinitionError, "<io.StringIO>:  1: External feature is implied above: in the following line:\nB external"):
            module.defs(io.StringIO("A implies B\nB external"))
        with self.assertRaisesRegex(module.DefinitionError, "<io.StringIO>:  0: <feature> implies \\[<feature>...\\] in the following line:\nA implies"):
            module.defs(io.StringIO("A implies "))
        with self.assertRaisesRegex(module.DefinitionError, "<io.StringIO>:  1: Implied feature B is already defined as external above: in the following line:\nA implies B,"):
            module.defs(io.StringIO("B external\nA implies B,"))
        with self.assertRaisesRegex(module.DefinitionError, "<io.StringIO>:  0: <feature> requires <expr> in the following line:\nA requires"):
            module.defs(io.StringIO("A requires "))

        # detect cyclic dependencies
        ruleset = "A implies B\nB implies A"
        with self.assertRaises(graphlib.CycleError):
            module.defs(io.StringIO(ruleset))
        ruleset = "A implies B\nB implies C\nC implies A"
        with self.assertRaises(graphlib.CycleError):
            module.defs(io.StringIO(ruleset))


if __name__ == "__main__":
    ut.main()
