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
# This module parses the feature definition file features.def
#
import fileinput, string, re

class SyntaxError:
    def __init__(self, message, instead):
        self.message = message
        self.filename = fileinput.filename()
        self.lineno = fileinput.filelineno()
        self.instead = instead
    def __str__(self):
        return '%s: %2d: %s in the following line:\n%s' % \
            (self.filename, self.lineno, self.message, self.instead)

def toCPPExpr(expr):
    expr = expr.replace('and', ' && ')
    expr = expr.replace('or', ' || ')
    expr = expr.replace('not', ' !')
    expr = re.sub('([A-Z0-9_]+)', 'defined(\\1)', expr)
    return expr

class defs:
    def __init__(self, filename):
        # complete set of all defined features
        features = set()
        # list of implications (pairs of feature -> implied feature)
        implications = list()
        # list of requirements (pairs of feature -> requirement expr)
        requirements = list()
        # list of derivations (pairs of feature -> derivation expr)
        derivations = list()
        # list of external features
        externals = set()

        for line in fileinput.input(filename):
            line = line.strip()
            # Ignore empty and comment lines
            if len(line) == 0 or line.startswith('#') \
                or line.startswith('//') or line.startswith('/*'): continue

            # Tokenify line
            tokens = line.split(None, 2)

            # Register the feature
            feature = tokens.pop(0)
            features.add(feature)

            # get the keyword
            if len(tokens) > 0:
                keyword = tokens.pop(0)
                if len(tokens) == 0:
                    rest = None
                else:
                    rest = tokens[0]

                # derived
                if keyword == 'equals':
                    if rest is None:
                        raise SyntaxError("<feature> equals <expr>", line)
                    derived = map((lambda(x,y,z):x), derivations)
                    if feature in derived:
                        raise SyntaxError("Derived feature is already defined above:", line);
                    if feature in externals:
                        raise SyntaxError("Derived feature is already defined as external above:", line);
                    features.discard(feature)
                    derivations.append((feature, rest, toCPPExpr(rest)))

                # externals
                elif keyword == 'external':
                    if rest is not None:
                        raise SyntaxError("<feature> external", line)
                    derived = map((lambda(x,y,z):x), derivations)
                    if feature in derived:
                        raise SyntaxError("External feature is already defined as derived above:", line);
                    implied = set(map((lambda (x,y):y), implications))
                    if feature in implied:
                        raise SyntaxError("External feature is implied above:", line);
                    externals.add(feature)
                    features.discard(feature)

                # implications
                elif keyword == 'implies':
                    if rest is None:
                        raise SyntaxError("<feature> implies [<feature>...]", line)
                    tokens = rest.split()
                    for implied in tokens:
                        if implied.endswith(','): implied = implied[:-1]
                        if implied in externals:
                            raise SyntaxError("Implied feature %s is already defined as external above:" % feature, line);

                        implications.append((feature, implied))

                # requires
                elif keyword == 'requires':
                    if rest is None:
                        raise SyntaxError("<feature> requires <expr>", line)
                    requirements.append((feature, rest, toCPPExpr(rest)))

        self.features = features
        self.requirements = requirements
        self.implications = implications
        self.derivations = derivations
        self.externals = externals

# Test whether all implied features or features in an expression are defined
