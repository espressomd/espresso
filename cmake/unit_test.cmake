#
# Copyright (C) 2016-2022 The ESPResSo project
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

# unit_test function
function(UNIT_TEST)
  cmake_parse_arguments(TEST "" "NAME;NUM_PROC" "SRC;DEPENDS" ${ARGN})
  add_executable(${TEST_NAME} ${TEST_SRC})
  # Build tests only when testing
  set_target_properties(${TEST_NAME} PROPERTIES EXCLUDE_FROM_ALL ON)
  target_link_libraries(${TEST_NAME} PRIVATE Boost::unit_test_framework)
  if(TEST_DEPENDS)
    target_link_libraries(${TEST_NAME} PRIVATE ${TEST_DEPENDS})
  endif()
  target_include_directories(${TEST_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src/core)
  target_link_libraries(${TEST_NAME} PRIVATE espresso::config espresso::cpp_flags)

  # If NUM_PROC is given, set up MPI parallel test case
  if(TEST_NUM_PROC)
    if(${TEST_NUM_PROC} GREATER ${ESPRESSO_TEST_NP})
      set(TEST_NUM_PROC ${ESPRESSO_TEST_NP})
    endif()
    espresso_set_mpiexec_tmpdir(${TEST_NAME})
    add_test(${TEST_NAME} ${MPIEXEC} ${ESPRESSO_MPIEXEC_OVERSUBSCRIBE}
             ${MPIEXEC_NUMPROC_FLAG} ${TEST_NUM_PROC} ${MPIEXEC_PREFLAGS}
             ${ESPRESSO_MPIEXEC_TMPDIR} ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}
             ${MPIEXEC_POSTFLAGS})
  else()
    add_test(${TEST_NAME} ${TEST_NAME})
  endif()

  if(ESPRESSO_WARNINGS_ARE_ERRORS)
    set(SANITIZERS_HALT_ON_ERROR "halt_on_error=1")
  else()
    set(SANITIZERS_HALT_ON_ERROR "halt_on_error=0")
  endif()
  list(APPEND TEST_ENV_VARIABLES "UBSAN_OPTIONS=suppressions=${CMAKE_SOURCE_DIR}/maintainer/CI/ubsan.supp:${SANITIZERS_HALT_ON_ERROR}:print_stacktrace=1")
  list(APPEND TEST_ENV_VARIABLES "ASAN_OPTIONS=${SANITIZERS_HALT_ON_ERROR}:detect_leaks=0:allocator_may_return_null=1")
  list(APPEND TEST_ENV_VARIABLES "MSAN_OPTIONS=${SANITIZERS_HALT_ON_ERROR}")
  if(NOT TEST_NUM_PROC AND ESPRESSO_MPIEXEC_GUARD_SINGLETON_NUMA AND "${TEST_DEPENDS}" MATCHES "(^|;)([Bb]oost::mpi|MPI::MPI_CXX)($|;)")
    list(APPEND TEST_ENV_VARIABLES "OMPI_MCA_hwloc_base_binding_policy=none")
  endif()
  set_tests_properties(
    ${TEST_NAME} PROPERTIES ENVIRONMENT "${TEST_ENV_VARIABLES}")

  add_dependencies(check_unit_tests ${TEST_NAME})
endfunction(UNIT_TEST)
