#
# Copyright (C) 2010-2022 The ESPResSo project
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

import sys
import subprocess

import defines
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
                # transposition
                d[(i, j)] = min(d[(i, j)], d[i - 2, j - 2] + cost)

    return d[lenstr1 - 1, lenstr2 - 1]


def handle_unknown(f, all_features):
    match = None
    max_dist = max(2, len(f) // 2)
    for d in all_features:
        dist = damerau_levenshtein_distance(f, d)
        if dist < max_dist:
            max_dist = dist
            match = d

    if match:
        return f"- unknown feature '{f}', did you mean '{match}'?"
    else:
        return f"- unknown feature '{f}'"


def check_myconfig(compiler, feature_file, myconfig, cmake_config=None):
    # query features from the compiler
    try:
        Defines = defines.Defines
        external_defs = []
        if cmake_config:
            external_features = Defines(compiler).defines(cmake_config)
            external_defs += ['-D' + s for s in external_features]
        my_features = Defines(compiler, flags=external_defs).defines(myconfig)
    except subprocess.CalledProcessError as ex:
        message = ex.output.decode("utf-8").split("\n")[0].strip()
        raise RuntimeError(
            f"Command `{' '.join(ex.cmd)}` returned non-zero exit code "
            f"{ex.returncode}, output: {message}.")

    # parse feature file
    defs = featuredefs.defs(feature_file)

    error_queue = []
    for e in (my_features & defs.externals):
        my_features.remove(e)
        error_queue.append(
            f"- external feature '{e}' cannot be defined in myconfig")

    for u in (my_features - defs.features):
        if u.startswith('__'):
            continue
        error_queue.append(handle_unknown(u, defs.features))

    if error_queue:
        error_report = "\n".join(error_queue)
        raise RuntimeError(
            f"There were errors in '{myconfig}':\n{error_report}")


if __name__ == "__main__":
    if len(sys.argv) > 4:
        cmake_config = sys.argv[4]
    else:
        cmake_config = None

    check_myconfig(sys.argv[1], sys.argv[2], sys.argv[3], cmake_config)
