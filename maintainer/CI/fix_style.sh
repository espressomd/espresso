#!/usr/bin/env sh
# Copyright (C) 2018-2019 The ESPResSo project
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

cd "$(git rev-parse --show-toplevel)"

CLANG_FORMAT_VER=6.0
AUTOPEP8_VER=1.3.4
PYCODESTYLE_VER=2.3.1

if ! git diff-index --quiet HEAD -- && [ "${1}" != "-f" ]; then
    echo "Warning, your working tree is not clean. Please commit your changes."
    echo "You can also call this script with the -f flag to proceed anyway, but"
    echo "you will then be unable to revert the formatting changes later."
    exit 1
fi

CLANGFORMAT="$(which clang-format-${CLANG_FORMAT_VER})"
if [ "${CLANGFORMAT}" = "" ]; then
    CLANGFORMAT="$(which clang-format)"
    if ! "${CLANGFORMAT}" --version | grep -qEo "version ${CLANG_FORMAT_VER}\.[0-9]+"; then
        echo "Could not find clang-format ${CLANG_FORMAT_VER}. ${CLANGFORMAT} is $(${CLANGFORMAT} --version | grep -Eo '[0-9\.]{5}' | head -n 1)."
        exit 2
    fi
fi

AUTOPEP8="$(which autopep8)"
if ! "${AUTOPEP8}" --version 2>&1 | grep -qFo "autopep8 ${AUTOPEP8_VER} (pycodestyle: ${PYCODESTYLE_VER})"; then
    echo "Could not find autopep8 ${AUTOPEP8_VER} with pycodestyle ${PYCODESTYLE_VER}"
    echo "${AUTOPEP8} is $(${AUTOPEP8} --version 2>&1)"
    exit 2
fi

find . \( -name '*.hpp' -o -name '*.cpp' -o -name '*.cu' -o -name '*.cuh' \) -not -path './libs/*' | xargs -r -n 5 -P 8 "${CLANGFORMAT}" -i -style=file || exit 3
find . \( -name '*.py' -o -name '*.pyx' -o -name '*.pxd' \) -not -path './libs/*' | xargs -r -n 5 -P 8 "${AUTOPEP8}" --ignore=E266,W291,W293 --in-place --aggressive || exit 3
find . -type f -executable ! -name '*.sh' ! -name '*.py' ! -name '*.sh.in' ! -name pypresso.cmakein -not -path './.git/*' | xargs -r -n 5 -P 8 chmod -x || exit 3

if [ "${CI}" != "" ]; then
    git --no-pager diff > style.patch
    maintainer/gh_post_style_patch.py || exit 1
fi
git diff-index --quiet HEAD --
if [ "${?}" = 1 ]; then
    if [ "${CI}" != "" ]; then
        echo "Failed style check. Download ${CI_JOB_URL}/artifacts/raw/style.patch to see which changes are necessary." >&2
    else
        echo "Failed style check" >&2
    fi
    exit 1
else
    echo "Passed style check"
fi

