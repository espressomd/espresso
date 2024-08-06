#
# Copyright (C) 2010-2022 The ESPResSo project
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

# - Find FFTW3
# Find the native FFTW3 includes and library, double precision
#
#  FFTW3_INCLUDE_DIR    - where to find fftw3.h
#  FFTW3_LIBRARIES   - List of libraries when using FFTW.
#  FFTW3_FOUND       - True if FFTW found.

if(FFTW3_INCLUDE_DIR)
  # Already in cache, be silent
  set(FFTW3_FIND_QUIETLY TRUE)
endif(FFTW3_INCLUDE_DIR)

find_path(FFTW3_INCLUDE_DIR fftw3.h)
find_library(FFTW3_LIBRARIES fftw3)
find_library(FFTW3F_LIBRARIES fftw3f)
find_path(FFTW3_MPI_INCLUDE_DIR fftw3-mpi.h)
find_library(FFTW3_MPI_LIBRARIES fftw3_mpi)
find_library(FFTW3F_MPI_LIBRARIES fftw3f_mpi)

# handle the QUIETLY and REQUIRED arguments and set FFTW3_FOUND to TRUE if all
# listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW3 DEFAULT_MSG FFTW3_LIBRARIES FFTW3F_LIBRARIES
                                  FFTW3_INCLUDE_DIR)
set(FPHSA_NAME_MISMATCHED 1)
find_package_handle_standard_args(FFTW3_MPI DEFAULT_MSG FFTW3_MPI_LIBRARIES FFTW3F_MPI_LIBRARIES
                                  FFTW3_MPI_INCLUDE_DIR)
unset(FPHSA_NAME_MISMATCHED)

mark_as_advanced(FFTW3_LIBRARIES FFTW3F_LIBRARIES FFTW3_INCLUDE_DIR FFTW3_MPI_LIBRARIES FFTW3F_MPI_LIBRARIES FFTW3_MPI_INCLUDE_DIR)


if(FFTW3_FOUND AND NOT TARGET FFTW3::FFTW3)
  add_library(FFTW3::FFTW3 INTERFACE IMPORTED)
  target_include_directories(FFTW3::FFTW3 INTERFACE "${FFTW3_INCLUDE_DIR}")
  target_link_libraries(FFTW3::FFTW3 INTERFACE "${FFTW3_LIBRARIES}" "${FFTW3F_LIBRARIES}")
endif()
