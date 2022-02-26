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
  target_link_libraries(${TEST_NAME} PRIVATE EspressoConfig cxx_interface)

  # If NUM_PROC is given, set up MPI parallel test case
  if(TEST_NUM_PROC)
    if(${TEST_NUM_PROC} GREATER ${TEST_NP})
      set(TEST_NUM_PROC ${TEST_NP})
    endif()
    set_mpiexec_tmpdir("${TEST_NAME}")
    add_test(${TEST_NAME} ${MPIEXEC} ${MPIEXEC_OVERSUBSCRIBE}
             ${MPIEXEC_NUMPROC_FLAG} ${TEST_NUM_PROC} ${MPIEXEC_PREFLAGS}
             ${MPIEXEC_TMPDIR} ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}
             ${MPIEXEC_POSTFLAGS})
  else()
    add_test(${TEST_NAME} ${TEST_NAME})
  endif()

  if(WARNINGS_ARE_ERRORS)
    set(SANITIZERS_HALT_ON_ERROR "halt_on_error=1")
  else()
    set(SANITIZERS_HALT_ON_ERROR "halt_on_error=0")
  endif()
  set(UBSAN_OPTIONS "UBSAN_OPTIONS=suppressions=${CMAKE_SOURCE_DIR}/maintainer/CI/ubsan.supp:${SANITIZERS_HALT_ON_ERROR}:print_stacktrace=1")
  set(ASAN_OPTIONS "ASAN_OPTIONS=${SANITIZERS_HALT_ON_ERROR}:detect_leaks=0:allocator_may_return_null=1")
  set(MSAN_OPTIONS "MSAN_OPTIONS=${SANITIZERS_HALT_ON_ERROR}")
  set_tests_properties(
    ${TEST_NAME} PROPERTIES ENVIRONMENT
                            "${UBSAN_OPTIONS} ${ASAN_OPTIONS} ${MSAN_OPTIONS}")

  add_dependencies(check_unit_tests ${TEST_NAME})
endfunction(UNIT_TEST)
