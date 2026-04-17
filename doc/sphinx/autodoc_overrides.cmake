#
# Copyright (C) 2023-2026 The ESPResSo project
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

file(READ ${FILEPATH} CONTENT)
set(TOKEN ".. automodule:: ${SUBMODULE}")
string(FIND "${CONTENT}" "${TOKEN}" POSITION)
if(POSITION EQUAL -1)
  message(FATAL_ERROR "Cannot find string ${TOKEN}")
endif()
string(REGEX REPLACE "${TOKEN}" "${TOKEN}\n   :member-order: bysource" CONTENT "${CONTENT}")
file(WRITE ${FILEPATH} "${CONTENT}")
