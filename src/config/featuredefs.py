#!/usr/bin/env python3
#
# Copyright (C) 2013-2026 The ESPResSo project
# Copyright (C) 2012 Olaf Lenz
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
This module parses the feature definition file :file:`features.def`.
"""
import fileinput
import graphlib
import io
import re


class DefinitionError(Exception):

    def __init__(self, message, context):
        self.message = message
        self.filename, self.lineno, self.instead = context

    def __str__(self):
        return f"{self.filename}: {self.lineno:>2}: {self.message} in the following line:\n{self.instead}"  # nopep8


def toCPPExpr(expr):
    expr = expr.replace("and", "&&")
    expr = expr.replace("or", "||")
    expr = expr.replace('not', "!")
    expr = re.sub("([A-Z0-9_]+)", r"defined(ESPRESSO_\g<1>)", expr)
    return expr


class defs:

    def __init__(self, filename):
        # complete set of all defined features
        allfeatures = set()
        # list of implications (pairs of feature -> implied feature)
        implications = list()
        # list of requirements (pairs of feature -> requirement expr)
        requirements = list()
        # set of derived features
        derived = set()
        # list of derivations (pairs of feature -> derivation expr)
        derivations = list()
        # list of external features
        externals = set()

        buffer = None
        if isinstance(filename, io.StringIO):
            buffer = filename
            filename = "<io.StringIO>"
        else:
            buffer = fileinput.input(filename)

        for lineno, line in enumerate(buffer):
            line = line.strip()
            context = (filename, lineno, line)

            # Ignore empty and comment lines
            if not line or line.startswith(('#', '//', '/*')):
                assert not line.startswith('/*') or line.endswith('*/')
                continue

            # Tokenify line
            tokens = line.split(None, 2)

            # Register the feature
            feature = tokens.pop(0)
            allfeatures.add(feature)

            # get the keyword
            if tokens:
                keyword = tokens.pop(0)
                if not tokens:
                    rest = None
                else:
                    rest = tokens[0]

                # derived
                if keyword == 'equals':
                    if rest is None:
                        raise DefinitionError(
                            "<feature> equals <expr>", context)
                    if feature in derived:
                        raise DefinitionError(
                            "Derived feature is already defined above:", context)
                    if feature in externals:
                        raise DefinitionError(
                            "Derived feature is already defined as external above:", context)
                    derived.add(feature)
                    derivations.append((feature, rest, toCPPExpr(rest)))

                # externals
                elif keyword == 'external':
                    if rest is not None:
                        raise DefinitionError("<feature> external", context)
                    if feature in derived:
                        raise DefinitionError(
                            "External feature is already defined as derived above:", context)
                    implied = set(map((lambda x_y: x_y[1]), implications))
                    if feature in implied:
                        raise DefinitionError(
                            "External feature is implied above:", context)
                    externals.add(feature)

                # implications
                elif keyword == 'implies':
                    if rest is None:
                        raise DefinitionError(
                            "<feature> implies [<feature>...]", context)
                    tokens = rest.split()
                    for implied in tokens:
                        if implied.endswith(','):
                            implied = implied[:-1]
                        if implied in externals:
                            raise DefinitionError(
                                f"Implied feature {implied} is already defined as external above:", context)

                        implications.append((feature, implied))

                # requires
                elif keyword == 'requires':
                    if rest is None:
                        raise DefinitionError(
                            "<feature> requires <expr>", context)
                    requirements.append((feature, rest, toCPPExpr(rest)))

        # allfeatures minus externals and derived
        features = allfeatures.difference(derived)
        features = features.difference(externals)

        self.allfeatures = allfeatures
        self.features = features
        self.requirements = requirements
        self.implications = implications
        self.derived = derived
        self.derivations = derivations
        self.externals = externals
        self.implications_topologically_ordered = self.process_implications()

    def process_implications(self):
        """
        Detect transitive chains in the graph of implication rules.
        Topological ordering is used to guarantee that any feature ``D`` is
        defined before any dependent rule ``D -> E`` is evaluated, such that
        a chain ``C -> D -> E`` works correctly when the preprocessor reads
        the feature config file from top to bottom. When multiple rules
        independently imply the same feature, they can be folded with syntax
        ``#if defined(A) or defined(B)`` to reduce the config file size.
        """
        graph = {}
        for feature, implied_feature in self.implications:
            graph.setdefault(feature, set())
            graph.setdefault(implied_feature, set()).add(feature)

        ts = graphlib.TopologicalSorter(graph)
        ts.prepare()
        rounds = []
        while ts.is_active():
            rules = {}
            # process rules that can be independently satisfied in this round
            for node in sorted(ts.get_ready()):
                for feature, implied_feature in self.implications:
                    if feature == node:
                        rules.setdefault(implied_feature, set()).add(feature)
                ts.done(node)
            if rules:
                rounds.append(rules)
        return rounds


class cmakedefs:
    def __init__(self, filename):
        self.externals = set()

        re_pattern = re.compile(
            "^#cmakedefine +ESPRESSO_BUILD_WITH_([A-Za-z0-9_]+) *$")
        inside_multiline_comment = False
        for line in fileinput.input(filename):
            line = line.strip()
            # Ignore empty and comment lines
            if inside_multiline_comment or line.startswith('/*'):
                inside_multiline_comment = not line.endswith('*/')
                continue
            if not line or line.startswith('//'):
                continue

            m = re_pattern.search(line)
            if m:
                self.externals.add(m.group(1))
