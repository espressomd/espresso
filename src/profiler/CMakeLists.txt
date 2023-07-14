#
# Copyright (C) 2019-2022 The ESPResSo project
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

add_library(espresso_profiler INTERFACE)
add_library(espresso::profiler ALIAS espresso_profiler)

if(ESPRESSO_BUILD_WITH_VALGRIND)
  target_include_directories(espresso_profiler
                             INTERFACE "${VALGRIND_INCLUDE_DIRS}")
endif()

if(ESPRESSO_BUILD_WITH_CALIPER)
  target_link_libraries(espresso_profiler INTERFACE caliper)
  # workaround to handle the case where the build directory is outside `_deps`
  cmake_path(GET CMAKE_BINARY_DIR FILENAME CMAKE_BINARY_DIR_BASENAME)
  set(caliper_BINARY_DIR_ALTERNATIVE
      "${CMAKE_BINARY_DIR}/${CMAKE_BINARY_DIR_BASENAME}/_deps/caliper-src")
  target_include_directories(
    espresso_profiler
    INTERFACE "${caliper_SOURCE_DIR}/include" "${caliper_BINARY_DIR}/include"
              "${caliper_BINARY_DIR_ALTERNATIVE}/include")
endif()
