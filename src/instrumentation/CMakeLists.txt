#
# Copyright (C) 2019-2026 The ESPResSo project
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

add_library(espresso_instrumentation INTERFACE)
add_library(espresso::instrumentation ALIAS espresso_instrumentation)
espresso_set_common_target_properties(espresso_instrumentation)
target_include_directories(
  espresso_instrumentation
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:include>)

if(ESPRESSO_BUILD_WITH_VALGRIND)
  target_include_directories(espresso_instrumentation
                             INTERFACE "${VALGRIND_INCLUDE_DIRS}")
endif()

if(ESPRESSO_BUILD_WITH_CALIPER)
  target_link_libraries(espresso_instrumentation INTERFACE caliper)
  target_include_directories(
    espresso_instrumentation INTERFACE "${caliper_SOURCE_DIR}/include"
                                       "${caliper_BINARY_DIR}/include")
endif()

if(ESPRESSO_BUILD_WITH_FPE)
  add_library(espresso_fpe SHARED src/fe_trap.cpp)
  add_library(espresso::fpe ALIAS espresso_fpe)
  espresso_set_common_target_properties(espresso_fpe)
  target_link_libraries(espresso_fpe PRIVATE espresso::config
                                             espresso::compiler_flags)
  target_include_directories(
    espresso_fpe PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                        $<INSTALL_INTERFACE:include>)
  target_link_libraries(espresso_instrumentation INTERFACE espresso::fpe)
  install(TARGETS espresso_fpe
          LIBRARY DESTINATION ${ESPRESSO_INSTALL_PYTHON}/espressomd)
endif()

add_subdirectory(tests)
