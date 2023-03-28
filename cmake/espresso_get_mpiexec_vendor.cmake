#
# Copyright (C) 2022 The ESPResSo project
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

# Find the MPI backend.
#
# This code sets the following variables:
#
#  ESPRESSO_MPIEXEC_VENDOR ESPRESSO_MPIEXEC_VERSION MPIEXEC_BACKEND_VERSION_REQUIRED

function(espresso_get_mpiexec_vendor)
  set(MPIEXEC_BACKEND_NAME "unknown")
  set(MPIEXEC_BACKEND_VERSION 0.0.0)
  set(MPIEXEC_BACKEND_VERSION_REQUIRED 0.0.0)
  execute_process(
    COMMAND ${MPIEXEC} --version RESULT_VARIABLE MPI_VERSION_RESULT
    OUTPUT_VARIABLE MPI_VERSION_OUTPUT ERROR_VARIABLE MPI_VERSION_OUTPUT)
  if(MPI_VERSION_RESULT EQUAL 0)
    if(MPI_VERSION_OUTPUT MATCHES "Intel\\(R\\) MPI Library")
      set(MPIEXEC_BACKEND_VERSION_REQUIRED 2021.0)
      set(MPIEXEC_BACKEND_NAME "Intel")
      string(REGEX REPLACE ".*Build ([0-9]+).*" "\\1" MPIEXEC_BACKEND_VERSION
                           ${MPI_VERSION_OUTPUT})
    endif()
    if(MPI_VERSION_OUTPUT MATCHES "HYDRA")
      set(MPIEXEC_BACKEND_VERSION_REQUIRED 3.4.1)
      set(MPIEXEC_BACKEND_NAME "MPICH")
      string(REGEX REPLACE ".*Version: +([0-9\\.]+).*" "\\1"
                           MPIEXEC_BACKEND_VERSION ${MPI_VERSION_OUTPUT})
    endif()
    if(MPI_VERSION_OUTPUT MATCHES "\\(Open(RTE| MPI)\\)")
      set(MPIEXEC_BACKEND_VERSION_REQUIRED 4.0.0)
      set(MPIEXEC_BACKEND_NAME "OpenMPI")
      string(REGEX REPLACE ".*\\(Open(RTE| MPI)\\) ([0-9\\.]+).*" "\\2"
                           MPIEXEC_BACKEND_VERSION ${MPI_VERSION_OUTPUT})
    endif()
  endif()
  set(ESPRESSO_MPIEXEC_VENDOR ${MPIEXEC_BACKEND_NAME} PARENT_SCOPE)
  set(ESPRESSO_MPIEXEC_VERSION ${MPIEXEC_BACKEND_VERSION} PARENT_SCOPE)
  set(ESPRESSO_MINIMAL_MPIEXEC_VERSION ${MPIEXEC_BACKEND_VERSION_REQUIRED} PARENT_SCOPE)
endfunction()
