#
# Copyright (C) 2016-2026 The ESPResSo project
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

# generate a test name based on the file name
function(espresso_unit_test_make_name)
  cmake_parse_arguments(TEST "" "SRC;NAME;SUFFIX" "" ${ARGN})
  if(NOT DEFINED TEST_NAME)
    cmake_path(GET TEST_SRC STEM TEST_NAME)
    set(TEST_NAME ${TEST_NAME} PARENT_SCOPE)
    if(DEFINED TEST_SUFFIX)
      set(TEST_NAME ${TEST_NAME}_${TEST_SUFFIX} PARENT_SCOPE)
    endif()
  elseif(DEFINED TEST_SUFFIX)
    message(FATAL_ERROR "Cannot provide both a NAME and SUFFIX (test: ${TEST_SRC})")
  endif()
endfunction()

# create the executable of a unit test
function(espresso_unit_test_executable)
  cmake_parse_arguments(TEST "" "SRC;NAME;SUFFIX;NUM_PROC;NUM_THREADS" "DEPENDS" ${ARGN})
  if(DEFINED TEST_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "Cannot parse arguments ${TEST_UNPARSED_ARGUMENTS}")
  endif()
  espresso_unit_test_make_name(${ARGV})
  add_executable(${TEST_NAME} ${TEST_SRC})
  espresso_set_common_target_properties(${TEST_NAME})
  # Build tests only when testing
  set_target_properties(${TEST_NAME} PROPERTIES EXCLUDE_FROM_ALL ON)
  target_link_libraries(
    ${TEST_NAME} PRIVATE Boost::unit_test_framework espresso::config
    espresso::compiler_flags espresso::tests::compiler_flags
    $<$<BOOL:${ESPRESSO_BUILD_WITH_FFTW}>:Heffte::Heffte>
    $<$<BOOL:${ESPRESSO_BUILD_WITH_WALBERLA}>:espresso::walberla>
    ${TEST_DEPENDS})
  target_include_directories(${TEST_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src/core)
  if(ESPRESSO_BUILD_WITH_CUDA)
    espresso_add_cuda_rpaths(${TEST_NAME}) # for GPU-aware MPI vendors
  endif()
endfunction()

# register the executable of a unit test in the CTest suite
function(espresso_unit_test_register)
  cmake_parse_arguments(TEST "" "SRC;NAME;SUFFIX;NUM_PROC;NUM_THREADS;TARGET" "DEPENDS" ${ARGN})
  if(DEFINED TEST_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "Cannot parse arguments ${TEST_UNPARSED_ARGUMENTS}")
  endif()
  espresso_unit_test_make_name(${ARGV})
  if(NOT DEFINED TEST_TARGET)
    set(TEST_TARGET ${TEST_NAME})
  endif()
  if(NOT DEFINED TEST_NUM_PROC)
    set(TEST_NUM_PROC 1)
  endif()
  if(TEST_NUM_PROC GREATER 1)
    if(${TEST_NUM_PROC} GREATER ${ESPRESSO_TEST_NP})
      set(TEST_NUM_PROC ${ESPRESSO_TEST_NP})
    endif()
    espresso_set_mpiexec_tmpdir(${TEST_NAME})
    add_test(NAME ${TEST_NAME} COMMAND ${MPIEXEC} ${MPIEXEC_PREFLAGS}
             ${ESPRESSO_MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG} ${TEST_NUM_PROC}
             ${ESPRESSO_MPIEXEC_TMPDIR} ${CMAKE_CURRENT_BINARY_DIR}/${TEST_TARGET}
             ${MPIEXEC_POSTFLAGS})
  else()
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_TARGET})
    if(ESPRESSO_MPIEXEC_GUARD_SINGLETON_NUMA)
      list(APPEND TEST_ENV_VARIABLES "OMPI_MCA_hwloc_base_binding_policy=none")
    endif()
  endif()
  if(NOT DEFINED TEST_NUM_THREADS)
    set(TEST_NUM_THREADS 2)
  endif()
  if(${TEST_NUM_THREADS} GREATER ${ESPRESSO_TEST_NT})
    set(TEST_NUM_THREADS ${ESPRESSO_TEST_NT})
  endif()
  list(APPEND TEST_ENV_VARIABLES "OMP_PROC_BIND=false" "OMP_NUM_THREADS=${TEST_NUM_THREADS}")

  if(ESPRESSO_WARNINGS_ARE_ERRORS)
    set(SANITIZERS_HALT_ON_ERROR "halt_on_error=1")
  else()
    set(SANITIZERS_HALT_ON_ERROR "halt_on_error=0")
  endif()
  list(APPEND TEST_ENV_VARIABLES "UBSAN_OPTIONS=suppressions=${CMAKE_SOURCE_DIR}/maintainer/CI/ubsan.supp:${SANITIZERS_HALT_ON_ERROR}:print_stacktrace=1")
  list(APPEND TEST_ENV_VARIABLES "ASAN_OPTIONS=${SANITIZERS_HALT_ON_ERROR}:detect_leaks=0:allocator_may_return_null=1")
  list(APPEND TEST_ENV_VARIABLES "MSAN_OPTIONS=${SANITIZERS_HALT_ON_ERROR}")
  math(EXPR TEST_NUM_CORES "${TEST_NUM_PROC} * ${TEST_NUM_THREADS}")
  set_tests_properties(
    ${TEST_NAME} PROPERTIES ENVIRONMENT "${TEST_ENV_VARIABLES}"
                            PROCESSORS ${TEST_NUM_CORES}
                            LABELS "unit_test")

  add_dependencies(unit_tests_executables ${TEST_TARGET})
endfunction()

# create the executable of a unit test and register it in the CTest suite
function(espresso_unit_test)
  espresso_unit_test_executable(${ARGV})
  espresso_unit_test_register(${ARGV})
endfunction()
