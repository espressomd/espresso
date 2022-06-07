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
#  MPIEXEC_BACKEND_NAME MPIEXEC_BACKEND_VERSION

set(MPIEXEC_BACKEND_NAME "unknown")
set(MPIEXEC_BACKEND_VERSION 0.0.0)

execute_process(
  COMMAND ${MPIEXEC} --version RESULT_VARIABLE mpi_version_result
  OUTPUT_VARIABLE mpi_version_output ERROR_VARIABLE mpi_version_output)
if(mpi_version_result EQUAL 0)
  if(mpi_version_output MATCHES "Intel\\(R\\) MPI Library")
    set(MPIEXEC_BACKEND_NAME "Intel")
    string(REGEX REPLACE ".*Build ([0-9]+).*" "\\1" MPIEXEC_BACKEND_VERSION ${mpi_version_output})
  endif()
  if(mpi_version_output MATCHES "HYDRA")
    set(MPIEXEC_BACKEND_NAME "MPICH")
    string(REGEX REPLACE ".*Version: +([0-9\\.]+).*" "\\1" MPIEXEC_BACKEND_VERSION ${mpi_version_output})
  endif()
  if(mpi_version_output MATCHES "\\(Open(RTE| MPI)\\)")
    set(MPIEXEC_BACKEND_NAME "OpenMPI")
    string(REGEX REPLACE ".*\\(Open(RTE| MPI)\\) ([0-9\\.]+).*" "\\2" MPIEXEC_BACKEND_VERSION ${mpi_version_output})
  endif()
endif()

include( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MpiexecBackend REQUIRED_VARS MPIEXEC)
