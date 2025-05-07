#
# Copyright (C) 2016-2025 The ESPResSo project
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

if(ESPRESSO_BUILD_WITH_HDF5)
  add_library(espresso_hdf5 SHARED h5md_core.cpp h5md_specification.cpp)
  add_library(espresso::hdf5 ALIAS espresso_hdf5)
  if(ESPRESSO_BUILD_WITH_CLANG_TIDY)
    set_target_properties(espresso_hdf5 PROPERTIES CXX_CLANG_TIDY
                                                   "${H5XX_CXX_CLANG_TIDY}")
  endif()
  target_include_directories(espresso_hdf5 PRIVATE ${CMAKE_SOURCE_DIR}/src/core)
  target_link_libraries(
    espresso_hdf5
    PRIVATE HighFive hdf5 espresso::cpp_flags MPI::MPI_CXX Boost::mpi
            Boost::filesystem espresso::utils espresso::config
            espresso::instrumentation)
  target_link_libraries(espresso_core PUBLIC espresso_hdf5)
  install(TARGETS espresso_hdf5
          LIBRARY DESTINATION ${ESPRESSO_INSTALL_PYTHON}/espressomd)
endif()
