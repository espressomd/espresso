#!/usr/bin/env bash
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

# list of commits to benchmark
commits="HEAD"

# list of directories to checkout
directories="../src ../libs"

cleanup() {
  # empty for now
  :
}

abort() {
  echo "An error occurred in suite.sh, exiting now" >&2
  echo "Command that failed: ${BASH_COMMAND}" >&2
  cleanup
  exit 1
}

trap abort EXIT
set -e

# move to top-level directory
cd "$(git rev-parse --show-toplevel)"
build_dir="$(realpath build-benchmarks)"

# move to build directory
if [ -d "${build_dir}" ]; then
  rm -rf "${build_dir}"
fi
mkdir "${build_dir}"
cd "${build_dir}"

# check for unstaged changes
if [ ! -z "$(git status --porcelain -- ${directories})" ]; then
  echo "fatal: you have unstaged changes, please commit or stash them:"
  git diff-index --name-only HEAD -- ${directories}
  exit 1
fi

cleanup() {
  # restore files in source directory
  git checkout HEAD -- ${directories}
}

# prepare output files
rm -f benchmarks.log
cat > benchmarks_suite.csv << EOF
"commit","config","script","arguments","cores","mean","ci","nsteps","duration"
EOF

# run benchmarks
for commit in ${commits}; do
  echo "### commit ${commit}" >> benchmarks.log
  git checkout ${commit} -- ${directories}
  bash ../maintainer/benchmarks/runner.sh
  sed -ri "s/^/\"${commit}\",/" benchmarks.csv
  tail -n +2 benchmarks.csv >> benchmarks_suite.csv
done

rm benchmarks.csv

trap : EXIT
cleanup
