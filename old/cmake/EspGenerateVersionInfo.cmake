# Copyright (C) 2009,2010 Christoph Junghans
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
macro(esp_generate_version_info ESPRESSO_VERSION ESPRESSO_TIMESTAMP) 
  execute_process( COMMAND sed -n "/^* /{s/.*\\(v[^[:space:]]*\\).*/\\1/p;q}" ${CMAKE_CURRENT_SOURCE_DIR}/RELEASE_NOTES OUTPUT_VARIABLE VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
  MESSAGE(STATUS "Found version ${VERSION}")
  SET(${ESPRESSO_VERSION} ${VERSION})
  execute_process( COMMAND sed -n "/^(/{s/^(...\\(.*\\)/\\1/p;q}" ${CMAKE_CURRENT_SOURCE_DIR}/RELEASE_NOTES OUTPUT_VARIABLE TIMESTAMP OUTPUT_STRIP_TRAILING_WHITESPACE)
  MESSAGE(STATUS "Last modified ${TIMESTAMP}")
  SET(${ESPRESSO_TIMESTAMP} ${TIMESTAMP})
endmacro(esp_generate_version_info ESPRESSO_VERSION ESPRESSO_TIMESTAMP) 
