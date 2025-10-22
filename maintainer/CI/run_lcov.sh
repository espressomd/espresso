#!/usr/bin/env sh
#
# Copyright (C) 2017-2024 The ESPResSo project
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

set -e

output="${1:-coverage.info}"
bindir="$(realpath .)"
srcdir="$(sed -nr "s/^ESPResSo_SOURCE_DIR:STATIC=(.+)/\1/p" "${bindir}/CMakeCache.txt")"

if [ "${srcdir}" = "" ]; then
  echo "Cannot extract ESPResSo_SOURCE_DIR variable from the CMake cache" >&2
  exit 2
fi

lcov --gcov-tool "${GCOV:-gcov}" \
     --quiet \
     --ignore-errors graph,mismatch,mismatch,gcov,unused \
     --directory . \
     --filter brace,blank,range,region \
     --capture \
     --rc lcov_json_module="JSON::XS" \
     --exclude "/usr/*" \
     --exclude "*/tmpxft_*cudafe1.stub.*" \
     --exclude "${bindir}/_deps/*" \
     --exclude "${bindir}/src/python/espressomd/*" \
     --exclude "${srcdir}/src/walberla_bridge/src/*/generated_kernels/*" \
     --exclude "${srcdir}/libs/*" \
     --output-file "${output}"
