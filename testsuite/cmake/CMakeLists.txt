if(NOT DEFINED TEST_NP)
  include(ProcessorCount)
  ProcessorCount(NP)
  math(EXPR TEST_NP "${NP}/2 + 1")
endif()

if(EXISTS ${MPIEXEC})
  # OpenMPI 3.0 and higher checks the number of processes against the number of CPUs
  execute_process(COMMAND ${MPIEXEC} --version RESULT_VARIABLE mpi_version_result OUTPUT_VARIABLE mpi_version_output ERROR_VARIABLE mpi_version_output)
  if (mpi_version_result EQUAL 0 AND mpi_version_output MATCHES "\\(Open(RTE| MPI)\\) ([3-9]\\.|1[0-9])")
    set(MPIEXEC_OVERSUBSCRIBE "-oversubscribe")
  else()
    set(MPIEXEC_OVERSUBSCRIBE "")
  endif()
endif()

function(CMAKE_TEST)
    cmake_parse_arguments(TEST "" "FILE;MAX_NUM_PROC;RUN_WITH_MPI" "DEPENDENCIES;CONFIGURATIONS" ${ARGN})
  get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
  configure_file(${TEST_FILE} ${CMAKE_CURRENT_BINARY_DIR}/${TEST_FILE} @ONLY)
  foreach(dependency IN LISTS TEST_DEPENDENCIES)
      configure_file(${dependency} ${CMAKE_CURRENT_BINARY_DIR}/${dependency} @ONLY)
  endforeach(dependency)
  set(TEST_FILE "${CMAKE_CURRENT_BINARY_DIR}/${TEST_FILE}")

  # Default values
  if (NOT DEFINED TEST_RUN_WITH_MPI)
    set(TEST_RUN_WITH_MPI TRUE)
  endif()
  set(TEST_CONFIGURATIONS ${TEST_CONFIGURATIONS} "serial")

  if(${TEST_MAX_NUM_PROC} LESS ${TEST_NP})
    set(TEST_NUM_PROC ${TEST_MAX_NUM_PROC})
  else()
    set(TEST_NUM_PROC ${TEST_NP})
  endif()

  if(EXISTS ${MPIEXEC} AND ${TEST_RUN_WITH_MPI})
    add_test(NAME ${TEST_NAME}
      COMMAND
      ${MPIEXEC} ${MPIEXEC_OVERSUBSCRIBE} ${MPIEXEC_NUMPROC_FLAG} ${TEST_NUM_PROC} ${TEST_FILE}
      CONFIGURATIONS
      ${TEST_CONFIGURATIONS}
      )
  else()
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_FILE})
  endif()

  set(cmake_tests ${cmake_tests} ${TEST_FILE} PARENT_SCOPE)
endfunction(CMAKE_TEST)

find_package(PythonInterp REQUIRED)
find_package(Doxygen)
find_package(Sphinx 1.6.6)
set(DESTDIR "${CMAKE_CURRENT_BINARY_DIR}/install")

cmake_test(FILE test_install.sh MAX_NUM_PROC 1 RUN_WITH_MPI FALSE DEPENDENCIES BashUnitTests.sh)
add_custom_target(setup_install COMMAND make install DESTDIR="${DESTDIR}" WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
cmake_test(FILE test_python_bindings.sh MAX_NUM_PROC 1 RUN_WITH_MPI FALSE DEPENDENCIES BashUnitTests.sh)

if(DOXYGEN_FOUND)
cmake_test(FILE test_doxygen.sh MAX_NUM_PROC 1 RUN_WITH_MPI FALSE DEPENDENCIES BashUnitTests.sh)
add_custom_target(setup_doxygen COMMAND make doc WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif(DOXYGEN_FOUND)

if(SPHINX_FOUND)
cmake_test(FILE test_sphinx.sh MAX_NUM_PROC 1 RUN_WITH_MPI FALSE DEPENDENCIES BashUnitTests.sh)
add_custom_target(setup_sphinx COMMAND make sphinx WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif(SPHINX_FOUND)

add_custom_target(check_cmake_serial COMMAND ${CMAKE_CTEST_COMMAND} -C serial --output-on-failure)
add_dependencies(check_cmake_serial setup_install)
if(DOXYGEN_FOUND)
add_dependencies(check_cmake_serial setup_doxygen)
endif(DOXYGEN_FOUND)
if(SPHINX_FOUND)
add_dependencies(check_cmake_serial setup_sphinx)
endif(SPHINX_FOUND)

add_custom_target(check_cmake)
add_dependencies(check_cmake check_cmake_serial)

add_dependencies(check check_cmake)
