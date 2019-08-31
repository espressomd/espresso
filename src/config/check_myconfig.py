# Copyright (C) 2010-2018 The ESPResSo project
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
import sys
from subprocess import CalledProcessError

from defines import Defines
import featuredefs


def damerau_levenshtein_distance(s1, s2):
    d = {}
    lenstr1 = len(s1)
    lenstr2 = len(s2)
    for i in range(-1, lenstr1 + 1):
        d[(i, -1)] = i + 1
    for j in range(-1, lenstr2 + 1):
        d[(-1, j)] = j + 1

    for i in range(lenstr1):
        for j in range(lenstr2):
            if s1[i] == s2[j]:
                cost = 0
            else:
                cost = 1
            d[(i, j)] = min(
                d[(i - 1, j)] + 1,  # deletion
                           d[(i, j - 1)] + 1,  # insertion
                           d[(i - 1, j - 1)] + cost,  # substitution
            )
            if i and j and s1[i] == s2[j - 1] and s1[i - 1] == s2[j]:
                d[(i, j)] = min(d[(i, j)], d[i - 2, j - 2] + cost)
                   # transposition

    return d[lenstr1 - 1, lenstr2 - 1]


def handle_unknown(f, all_features):
    match = None
    max_dist = max(2, len(f) // 2)
    for d in all_features:
        dist = damerau_levenshtein_distance(f, d)
        if dist < max_dist:
            min_dist = dist
            match = d

    if match:
        print("Unknown feature '{}', did you mean '{}'?".format(f, match))
    else:
        print("Unknown feature '{}'".format(f))


class FeatureError(Exception):
    pass


def print_exception(ex):
    print("""Skipped external header because {} returned non-zero exit code {},
             output: {}.""".format(' '.join(ex.cmd), ex.returncode, ex.output.strip()))


def check_myconfig(compiler, feature_file, myconfig, pre_header=None):
    # This does not work on all compilers, so if the parsing fails
    # we just bail out.
    external_defs = []

    if pre_header:
        try:
            external_features = Defines(compiler).defines(pre_header)
        except CalledProcessError as ex:
            print_exception(ex)
            return

        external_defs = ['-D' + s for s in external_features]

    try:
        my_features = Defines(compiler, flags=external_defs).defines(myconfig)
    except CalledProcessError as ex:
        print_exception(ex)
        return

    # Parse feature file
    defs = featuredefs.defs(feature_file)

    error_state = False
    for e in (my_features & defs.externals):
        error_state = True
        my_features.remove(e)
        print(
            "External feature '{}' can not be defined in myconfig.".format(e))

    for u in (my_features - defs.features):
        if u.startswith('__'):
            continue
        error_state = True
        handle_unknown(u, defs.features)

    if error_state:
        raise FeatureError("There were errors in '{}'".format(sys.argv[3]))
    else:
        return

if __name__ == "__main__":
    if len(sys.argv) > 4:
        pre_header = sys.argv[4]
    else:
        pre_header = None

    try:
        check_myconfig(sys.argv[1], sys.argv[2], sys.argv[3], pre_header)
        sys.exit()
    except FeatureError:
        sys.exit("There were errors in '{}'".format(sys.argv[3]))
