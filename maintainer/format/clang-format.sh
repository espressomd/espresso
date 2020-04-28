#!/bin/sh
# Copyright (C) 2018-2020 The ESPResSo project
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

CLANG_FORMAT_VER=9.0
if hash clang-format-${CLANG_FORMAT_VER} 2>/dev/null; then
    CLANGFORMAT="$(which clang-format-${CLANG_FORMAT_VER})"
elif hash clang-format-${CLANG_FORMAT_VER%.*} 2>/dev/null; then
    CLANGFORMAT="$(which clang-format-${CLANG_FORMAT_VER%.*})"
elif hash clang-format 2>/dev/null; then
    CLANGFORMAT="$(which clang-format)"
else
    echo "No clang-format found."
    exit 2
fi

if ! "${CLANGFORMAT}" --version | grep -qEo "version ${CLANG_FORMAT_VER}\.[0-9]+"; then
    echo "Could not find clang-format ${CLANG_FORMAT_VER}. ${CLANGFORMAT} is $(${CLANGFORMAT} --version | grep -Eo '[0-9\.]{5}' | head -n 1)."
    exit 2
fi

${CLANGFORMAT} "$@"
