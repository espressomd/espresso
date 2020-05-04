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

CMAKE_FORMAT_VER=0.6.9
python3 -m cmake_format 2>&1 > /dev/null
if [ "$?" = "0" ]; then
    CMAKE_FORMAT="python3 -m cmake_format"
else
    echo "No cmake-format found."
    exit 2
fi


if ! ${CMAKE_FORMAT} --version | grep -qEo "${CMAKE_FORMAT_VER}"; then
    echo "Could not find cmake-format ${CMAKE_FORMAT_VER}."
    exit 2
fi

${CMAKE_FORMAT} "$@"
