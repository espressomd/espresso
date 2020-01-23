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


# Add copyright header to C++, CUDA, Doxygen, Python and shell files.

cd "$(git rev-parse --show-toplevel)"

# Get files without headers
all_files=$(maintainer/files_with_header.sh)
files_to_check=$(grep -iL copyright ${all_files})
py_files=$(echo ${files_to_check} | tr " " "\n" | grep -P '\.(pyx?|pxd|sh)$')
cpp_files=$(echo ${files_to_check} | tr " " "\n" | grep -P '\.([c|h]pp|cuh?|dox)$')

tmp=$(mktemp)
# process python/cython/bash files
for f in ${py_files}; do
  head -n1 "${f}" | grep -q '^#!'
  if [ "${?}" = 0 ]; then
    # preserve shebang on first line
    (head -n1 "${f}"
     sed -e 's/^/# /' maintainer/header_template.txt | sed 's/ $//'
     tail -n+2 "${f}") > "${tmp}"
  else
    (sed -e 's/^/# /' maintainer/header_template.txt | sed 's/ $//'; cat "${f}") > "${tmp}"
  fi
  cp "${tmp}" "${f}"
  echo "${f}"
done
# process c++/cuda/doxygen files
for f in ${cpp_files}; do
  (echo '/*'
   sed -e 's/^/ * /' maintainer/header_template.txt | sed 's/ $//'
   echo ' */'
   cat "${f}") > "${tmp}"
  cp "${tmp}" "${f}"
  echo "${f}"
done
rm "${tmp}"
