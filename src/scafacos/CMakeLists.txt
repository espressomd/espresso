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

add_library(Espresso_scafacos SHARED src/Scafacos.cpp src/Coulomb.cpp
                                     src/Dipoles.cpp)
add_library(Espresso::scafacos ALIAS Espresso_scafacos)
if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 10 AND INSIDE_DOCKER
   AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS_EQUAL 11 AND WITH_COVERAGE)
  target_link_libraries(Espresso_scafacos
                        PRIVATE "-L/usr/lib/gcc/x86_64-linux-gnu/10")
endif()
target_link_libraries(Espresso_scafacos PUBLIC MPI::MPI_CXX
                      PRIVATE ${SCAFACOS_LDFLAGS} Espresso::cpp_flags)

target_include_directories(
  Espresso_scafacos
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
         $<INSTALL_INTERFACE:include>)
target_include_directories(Espresso_scafacos SYSTEM
                           PUBLIC ${SCAFACOS_INCLUDE_DIRS})

install(TARGETS Espresso_scafacos DESTINATION ${PYTHON_INSTDIR}/espressomd)
