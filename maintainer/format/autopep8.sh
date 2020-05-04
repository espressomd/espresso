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


AUTOPEP8_VER=1.3.4
PYCODESTYLE_VER=2.3.1

python3 -m autopep8 --help 2>&1 > /dev/null
if [ "$?" = "0" ]; then
    AUTOPEP8="python3 -m autopep8"
else
    echo "No autopep8 found."
    exit 2
fi

if ! ${AUTOPEP8} --version 2>&1 | grep -qFo "autopep8 ${AUTOPEP8_VER} (pycodestyle: ${PYCODESTYLE_VER})"; then
    echo "Could not find autopep8 ${AUTOPEP8_VER} with pycodestyle ${PYCODESTYLE_VER}"
    echo "${AUTOPEP8} is $(${AUTOPEP8} --version 2>&1)"
    exit 2
fi

${AUTOPEP8} "$@"
