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

add_library(espresso_utils INTERFACE)
add_library(espresso::utils ALIAS espresso_utils)
add_library(espresso_utils_mpi INTERFACE)
add_library(espresso::utils::mpi ALIAS espresso_utils)
target_include_directories(
  espresso_utils
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(espresso_utils INTERFACE Boost::serialization)
target_link_libraries(espresso_utils_mpi INTERFACE espresso::utils Boost::mpi
                                                   MPI::MPI_CXX)

if(ESPRESSO_BUILD_TESTS)
  add_subdirectory(tests)
endif()
