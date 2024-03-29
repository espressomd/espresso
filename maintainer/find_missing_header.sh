#!/usr/bin/env sh
#
# Copyright (C) 2012-2022 The ESPResSo project
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


# Check for missing GPL and copyright headers

files=$(maintainer/files_with_header.sh)
num_files=$(echo ${files} | wc -w)
current_year=$(date +%Y)

echo "Examining ${num_files} files."
echo

missing_current_copyright=$(grep -ILE "Copyright.*${current_year}" ${files})

echo "Copyright disclaimer missing the current year ${current_year}"
echo "--------------------------------------------------"
if [ -n "${missing_current_copyright}" ]; then
    grep -IlE "Copyright" ${missing_current_copyright}
else
    echo "NONE"
fi
echo

echo "Missing copyright disclaimer"
echo "----------------------------"
if [ -n "${missing_current_copyright}" ]; then
    grep -ILE "Copyright" ${missing_current_copyright}
else
    echo "NONE"
fi
echo

echo "Missing GPL/simple header"
echo "-------------------------"
nolicense=$(grep -ILE "((ESPResSo|This program) is free software|Copying and distribution of this file)" ${files})
if [ -n "${nolicense}" ]; then
    ls -1 ${nolicense}
else
    echo "NONE"
fi
echo

if [ -n "${nolicense}" ] || [ -n "${missing_current_copyright}" ]; then
  exit 1
fi
