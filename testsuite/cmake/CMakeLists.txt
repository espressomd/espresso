function(CMAKE_TEST)
    cmake_parse_arguments(TEST "" "FILE" "DEPENDENCIES" ${ARGN})
  get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
  configure_file(${TEST_FILE} ${CMAKE_CURRENT_BINARY_DIR}/${TEST_FILE} @ONLY)
  foreach(dependency IN LISTS TEST_DEPENDENCIES)
      configure_file(${dependency} ${CMAKE_CURRENT_BINARY_DIR}/${dependency} @ONLY)
  endforeach(dependency)
  set(TEST_FILE "${CMAKE_CURRENT_BINARY_DIR}/${TEST_FILE}")
  add_test(NAME ${TEST_NAME} COMMAND ${TEST_FILE})
  set_tests_properties(${TEST_NAME} PROPERTIES RUN_SERIAL ON)
  set(cmake_tests ${cmake_tests} ${TEST_FILE} PARENT_SCOPE)
endfunction(CMAKE_TEST)

execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c "from distutils import sysconfig as sc;print(sc.get_python_lib(prefix='${CMAKE_INSTALL_PREFIX}', plat_specific=True, standard_lib=False))"
  OUTPUT_VARIABLE PYTHON_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE)

cmake_test(FILE test_install.sh DEPENDENCIES BashUnitTests.sh)

add_custom_target(setup_install COMMAND make install WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_custom_target(check_cmake_install COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${TEST_TIMEOUT} -C serial --output-on-failure)
add_dependencies(check_cmake_install setup_install)

