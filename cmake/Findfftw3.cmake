#
# Copyright (C) 2025 The ESPResSo project
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

# Find the fftw3 libraries (with variants MPI, OpenMP, etc.) and header files.
#
# The following variables and imported targets are created:
#
# - `fftw3_INCLUDE_DIR`     : folder containing `fftw3.h`
# - `fftw3_MPI_INCLUDE_DIR` : folder containing `fftw3-mpi.h`
# - `fftw3_FOUND`           : true if `espresso::fftw3` exists and all required components were found
# - `fftw3_mpi_FOUND`       : true if `espresso::fftw3_mpi` exists
# - `fftw3_omp_FOUND`       : true if `espresso::fftw3_omp` exists
# - `fftw3_threads_FOUND`   : true if `espresso::fftw3_threads` exists

include(FindPackageHandleStandardArgs)

set(fftw3_PRECISIONS "" "f")
set(fftw3_COMPONENTS "mpi" "omp" "threads")
set(fftw3_MODULES ${fftw3_COMPONENTS})
list(TRANSFORM fftw3_MODULES PREPEND "_")
list(INSERT fftw3_MODULES 0 "")

find_path(fftw3_INCLUDE_DIR fftw3.h)
find_path(fftw3_MPI_INCLUDE_DIR fftw3-mpi.h)

foreach(PRECISION IN LISTS fftw3_PRECISIONS)
  foreach(COMPONENT IN LISTS fftw3_COMPONENTS)
    find_library(fftw3${PRECISION}_${COMPONENT}_LIBRARIES
                 fftw3${PRECISION}_${COMPONENT})
    mark_as_advanced(fftw3${PRECISION}_${COMPONENT}_LIBRARIES)
  endforeach()
  find_library(fftw3${PRECISION}_LIBRARIES fftw3${PRECISION})
  mark_as_advanced(fftw3${PRECISION}_LIBRARIES)
endforeach()

# dependencies bookkeeping
set(fftw3_LIBRARIES_REQUIRED fftw3_INCLUDE_DIR)
foreach(PRECISION IN LISTS fftw3_PRECISIONS)
  list(APPEND fftw3_LIBRARIES_REQUIRED fftw3${PRECISION}_LIBRARIES)
endforeach()

# find library components
set(FPHSA_NAME_MISMATCHED 1)
foreach(COMPONENT IN LISTS fftw3_COMPONENTS)
  if(${COMPONENT} IN_LIST fftw3_FIND_COMPONENTS)
    set(fftw3_COMPONENT_DEPS fftw3_INCLUDE_DIR)
    if(${COMPONENT} STREQUAL "mpi")
      list(APPEND fftw3_COMPONENT_DEPS fftw3_MPI_INCLUDE_DIR)
    endif()
    foreach(PRECISION IN LISTS fftw3_PRECISIONS)
      list(APPEND fftw3_COMPONENT_DEPS fftw3${PRECISION}_${COMPONENT}_LIBRARIES)
    endforeach()
    find_package_handle_standard_args(fftw3_${COMPONENT} DEFAULT_MSG ${fftw3_COMPONENT_DEPS})
    if(${fftw3_FIND_REQUIRED_${COMPONENT}})
      list(APPEND fftw3_LIBRARIES_REQUIRED ${fftw3_COMPONENT_DEPS})
    endif()
  endif()
endforeach()
list(REMOVE_DUPLICATES fftw3_LIBRARIES_REQUIRED)
unset(FPHSA_NAME_MISMATCHED)

# find main library, but fail if required components are missing
find_package_handle_standard_args(fftw3 DEFAULT_MSG ${fftw3_LIBRARIES_REQUIRED})

foreach(MODULE IN LISTS fftw3_MODULES)
  if(fftw3${MODULE}_FOUND AND NOT TARGET espresso::fftw3${MODULE})
    add_library(espresso::fftw3${MODULE} INTERFACE IMPORTED)
    target_include_directories(
      espresso::fftw3${MODULE}
      INTERFACE "${fftw3_INCLUDE_DIR}" "$<$<STREQUAL:${MODULE},_mpi>:${fftw3_MPI_INCLUDE_DIR}>")
    foreach(PRECISION IN LISTS fftw3_PRECISIONS)
      target_link_libraries(
        espresso::fftw3${MODULE}
        INTERFACE "${fftw3${PRECISION}${MODULE}_LIBRARIES}")
    endforeach()
  endif()
endforeach()
