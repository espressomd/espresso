#!/usr/bin/env sh
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

cd "$(git rev-parse --show-toplevel)" || exit 1

if ! git diff-index --quiet HEAD -- && [ "${1}" != "-f" ]; then
    echo "Warning, your working tree is not clean. Please commit your changes."
    echo "You can also call this script with the -f flag to proceed anyway, but"
    echo "you will then be unable to revert the formatting changes later."
    exit 1
fi

maintainer/lint/pre_commit.sh run --all-files
pre_commit_return_code="${?}"

git diff-index --quiet HEAD --
git_diff_return_code="${?}"

pylint_errors=0
if [ -f "pylint.log" ]; then
    pylint_errors=$(grep -Pc '^[a-z]+/.+?.py:[0-9]+:[0-9]+: [CRWEF][0-9]+:' "pylint.log")
fi

if [ -n "${CI}" ]; then
    git --no-pager diff > style.patch
    maintainer/gh_post_style_patch.py || exit 1
    maintainer/gh_post_pylint.py "${pylint_errors}" pylint.log || exit 1
    if [ "${git_diff_return_code}" != 0 ]; then
        echo "Failed style check. Download ${CI_JOB_URL}/artifacts/raw/style.patch to see which changes are necessary." >&2
    fi
    if [ "${pylint_errors}" != 0 ]; then
        echo "Failed pylint check: ${pylint_errors} errors" >&2
    fi
else
    if [ "${git_diff_return_code}" != 0 ]; then
        echo "Failed style check." >&2
    fi
    if [ "${pylint_errors}" != 0 ]; then
        echo "Failed pylint check." >&2
    fi
fi


exit ${pre_commit_return_code}
