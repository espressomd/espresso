function(BENCHMARK_TEST)
  python_scripts_test(${ARGV} TYPE benchmark) # forward arguments
  set(benchmarks_tests ${benchmarks_tests} ${TEST_FILE_CONFIGURED} PARENT_SCOPE)
endfunction(BENCHMARK_TEST)

# configure Python module importlib_wrapper.py
add_test(importlib_wrapper ${CMAKE_BINARY_DIR}/pypresso
         ${TEST_FILE_CONFIGURED_IMPORTLIB_WRAPPER})
set_tests_properties(importlib_wrapper PROPERTIES FIXTURES_SETUP
                                                  IMPORTLIB_WRAPPER)
set(benchmarks_tests ${benchmarks_tests}
                     ${TEST_FILE_CONFIGURED_IMPORTLIB_WRAPPER} PARENT_SCOPE)
configure_file(../importlib_wrapper.py
               ${CMAKE_CURRENT_BINARY_DIR}/importlib_wrapper.py)

# create a local copy of the benchmarks
set(BENCHMARKS_DIR "${CMAKE_CURRENT_BINARY_DIR}/local_benchmarks")
add_custom_target(
  local_benchmarks COMMAND ${CMAKE_COMMAND} -E remove_directory
                           ${BENCHMARKS_DIR} # cleanup
  COMMAND ${CMAKE_COMMAND} -E copy_directory
          ${CMAKE_SOURCE_DIR}/maintainer/benchmarks ${BENCHMARKS_DIR})

benchmark_test(FILE test_lj.py)
benchmark_test(FILE test_p3m.py)

add_custom_target(
  check_benchmarks COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${TEST_TIMEOUT}
                           -j${TEST_NP} $(ARGS) --output-on-failure)

add_dependencies(check_benchmarks pypresso local_benchmarks)
