include(ProcessorCount)
ProcessorCount(NP)

if(EXISTS ${MPIEXEC})
  # OpenMPI 3.0 and higher checks the number of processes against the number of
  # CPUs
  execute_process(
    COMMAND ${MPIEXEC} --version RESULT_VARIABLE mpi_version_result
    OUTPUT_VARIABLE mpi_version_output ERROR_VARIABLE mpi_version_output)
  if(mpi_version_result EQUAL 0 AND mpi_version_output MATCHES
                                    "\\(Open(RTE| MPI)\\) ([3-9]\\.|1[0-9])")
    set(MPIEXEC_OVERSUBSCRIBE "-oversubscribe")
  else()
    set(MPIEXEC_OVERSUBSCRIBE "")
  endif()
endif()

function(PYTHON_BENCHMARK)
  cmake_parse_arguments(
    BENCHMARK "" "FILE;RUN_WITH_MPI;MIN_NUM_PROC;MAX_NUM_PROC"
    "ARGUMENTS;DEPENDENCIES" ${ARGN})
  get_filename_component(BENCHMARK_NAME ${BENCHMARK_FILE} NAME_WE)
  foreach(argument IN LISTS BENCHMARK_ARGUMENTS)
    string(REGEX REPLACE "[^-a-zA-Z0-9_\\.]+" "_" argument ${argument})
    string(REGEX REPLACE "^[-_]+" "" argument ${argument})
    set(BENCHMARK_NAME "${BENCHMARK_NAME}__${argument}")
  endforeach(argument)
  configure_file(${BENCHMARK_FILE}
                 ${CMAKE_CURRENT_BINARY_DIR}/${BENCHMARK_FILE})
  foreach(dependency IN LISTS BENCHMARK_DEPENDENCIES)
    configure_file(${dependency} ${CMAKE_CURRENT_BINARY_DIR}/${dependency})
  endforeach(dependency)
  set(BENCHMARK_FILE "${CMAKE_CURRENT_BINARY_DIR}/${BENCHMARK_FILE}")
  list(APPEND BENCHMARK_ARGUMENTS
       "--output=${CMAKE_BINARY_DIR}/benchmarks.csv.part")

  # default values
  if(NOT DEFINED BENCHMARK_RUN_WITH_MPI)
    set(BENCHMARK_RUN_WITH_MPI TRUE)
  endif()
  if(NOT DEFINED BENCHMARK_MIN_NUM_PROC)
    set(BENCHMARK_MIN_NUM_PROC 1)
  endif()
  if(NOT DEFINED BENCHMARK_MAX_NUM_PROC)
    set(BENCHMARK_MAX_NUM_PROC ${NP})
  endif()
  # parallel schemes
  if(EXISTS ${MPIEXEC} AND ${BENCHMARK_RUN_WITH_MPI})
    set(BENCHMARK_CONFIGURATIONS "sentinel")
    foreach(nproc 1 2 4 8 16)
      if(${BENCHMARK_MAX_NUM_PROC} GREATER_EQUAL ${nproc}
         AND ${BENCHMARK_MIN_NUM_PROC} LESS_EQUAL ${nproc}
         AND ${NP} GREATER_EQUAL ${nproc})
        list(APPEND BENCHMARK_CONFIGURATIONS ${nproc})
      endif()
    endforeach(nproc)
    list(REMOVE_AT BENCHMARK_CONFIGURATIONS 0)
    foreach(nproc IN LISTS BENCHMARK_CONFIGURATIONS)
      set(BENCHMARK_TEST_NAME benchmark__${BENCHMARK_NAME}__parallel_${nproc})
      add_test(
        NAME ${BENCHMARK_TEST_NAME}
        COMMAND
          ${MPIEXEC} ${MPIEXEC_OVERSUBSCRIBE} ${MPIEXEC_NUMPROC_FLAG} ${nproc}
          ${MPIEXEC_PREFLAGS} ${CMAKE_BINARY_DIR}/pypresso ${BENCHMARK_FILE}
          ${BENCHMARK_ARGUMENTS} ${MPIEXEC_POSTFLAGS})
      set_tests_properties(${BENCHMARK_TEST_NAME} PROPERTIES RUN_SERIAL TRUE)
    endforeach(nproc)
  else()
    set(BENCHMARK_TEST_NAME benchmark__${BENCHMARK_NAME}__serial)
    add_test(NAME ${BENCHMARK_TEST_NAME}
             COMMAND ${CMAKE_BINARY_DIR}/pypresso ${BENCHMARK_FILE}
                     ${BENCHMARK_ARGUMENTS})
    set_tests_properties(${BENCHMARK_TEST_NAME} PROPERTIES RUN_SERIAL TRUE)
  endif()
endfunction(PYTHON_BENCHMARK)

python_benchmark(FILE lb.py ARGUMENTS "--box_l=32;--single_precision")
python_benchmark(FILE lb.py ARGUMENTS "--box_l=32")
python_benchmark(FILE lb.py ARGUMENTS "--box_l=64;--single_precision")
python_benchmark(FILE lb.py ARGUMENTS "--box_l=64")
python_benchmark(FILE lb.py ARGUMENTS "--box_l=128;--single_precision")
python_benchmark(FILE lb.py ARGUMENTS "--box_l=128")
if(NOT WALBERLA_USE_AVX)
  python_benchmark(FILE lb.py ARGUMENTS "--box_l=196;--single_precision")
  python_benchmark(FILE lb.py ARGUMENTS "--box_l=196")
endif()

add_custom_target(
  benchmarks_data
  COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks.py
          ${CMAKE_CURRENT_BINARY_DIR})
add_custom_target(
  benchmark_python COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${TEST_TIMEOUT}
                           ${CTEST_ARGS} --output-on-failure)

add_dependencies(benchmark benchmark_python benchmarks_data)
