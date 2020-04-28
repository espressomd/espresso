#!/usr/bin/env sh
# Copyright (C) 2017-2019 The ESPResSo project
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


# Detect python tests that are not included in CMakeLists.txt. Run this script
# in the testsuite folders for python, sample, tutorial and benchmark tests.

for T in *.py; do
    if grep -qF " ${T}" CMakeLists.txt; then
        continue
    else
        tracked_status=$(git ls-files -- "${T}")
        if [ -z "${tracked_status}" ]; then
            echo "File '${T}' is not tracked."
            continue
        else
            echo "File '${T}' is missing in CMakeLists.txt."
        fi
    fi
done
