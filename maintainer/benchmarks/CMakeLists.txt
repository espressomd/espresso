if(NOT DEFINED TEST_NP)
  include(ProcessorCount)
  processorcount(NP)
  math(EXPR TEST_NP "${NP}/2 + 1")
endif()

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
    set(BENCHMARK_CONFIGURATIONS "0")
    if(${NP} GREATER 0 AND ${BENCHMARK_MAX_NUM_PROC} GREATER 0
       AND ${BENCHMARK_MIN_NUM_PROC} LESS 2)
      list(APPEND BENCHMARK_CONFIGURATIONS 1)
    endif()
    if(${NP} GREATER 1 AND ${BENCHMARK_MAX_NUM_PROC} GREATER 1
       AND ${BENCHMARK_MIN_NUM_PROC} LESS 3)
      list(APPEND BENCHMARK_CONFIGURATIONS 2)
    endif()
    if(${NP} GREATER 3 AND ${BENCHMARK_MAX_NUM_PROC} GREATER 3
       AND ${BENCHMARK_MIN_NUM_PROC} LESS 5)
      list(APPEND BENCHMARK_CONFIGURATIONS 4)
    endif()
    if(${NP} GREATER 7 AND ${BENCHMARK_MAX_NUM_PROC} GREATER 7
       AND ${BENCHMARK_MIN_NUM_PROC} LESS 9)
      list(APPEND BENCHMARK_CONFIGURATIONS 8)
    endif()
    if(${NP} GREATER 15 AND ${BENCHMARK_MAX_NUM_PROC} GREATER 15
       AND ${BENCHMARK_MIN_NUM_PROC} LESS 17)
      list(APPEND BENCHMARK_CONFIGURATIONS 16)
    endif()
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

python_benchmark(FILE lj.py ARGUMENTS
                 "--particles_per_core=1000;--volume_fraction=0.50")
python_benchmark(FILE lj.py ARGUMENTS
                 "--particles_per_core=1000;--volume_fraction=0.02")
python_benchmark(FILE lj.py ARGUMENTS
                 "--particles_per_core=10000;--volume_fraction=0.50")
python_benchmark(FILE lj.py ARGUMENTS
                 "--particles_per_core=10000;--volume_fraction=0.02")
python_benchmark(
  FILE lj.py ARGUMENTS
  "--particles_per_core=1000;--volume_fraction=0.10;--bonds" RUN_WITH_MPI FALSE)
python_benchmark(
  FILE lj.py ARGUMENTS
  "--particles_per_core=10000;--volume_fraction=0.10;--bonds" RUN_WITH_MPI
  FALSE)
python_benchmark(
  FILE p3m.py ARGUMENTS
  "--particles_per_core=1000;--volume_fraction=0.25;--prefactor=4")
python_benchmark(
  FILE p3m.py ARGUMENTS
  "--particles_per_core=10000;--volume_fraction=0.25;--prefactor=4")

add_custom_target(
  benchmark_python COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${TEST_TIMEOUT}
                           $(ARGS) --output-on-failure)

add_dependencies(benchmark benchmark_python)
