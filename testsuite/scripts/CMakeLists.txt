set(TEST_FILE_CONFIGURED_IMPORTLIB_WRAPPER
    ${CMAKE_CURRENT_BINARY_DIR}/test_importlib_wrapper.py)
configure_file(importlib_wrapper.py
               ${CMAKE_CURRENT_BINARY_DIR}/importlib_wrapper.py)
configure_file(test_importlib_wrapper.py
               ${TEST_FILE_CONFIGURED_IMPORTLIB_WRAPPER})

macro(PYTHON_SCRIPTS_TEST)
  cmake_parse_arguments(TEST "" "FILE;SUFFIX;TYPE" "DEPENDENCIES;LABELS"
                        ${ARGN})
  get_filename_component(TEST_NAME ${TEST_FILE} NAME_WE)
  if(TEST_SUFFIX)
    set(TEST_NAME "${TEST_NAME}_with_${TEST_SUFFIX}")
  endif()
  set(TEST_FILE_CONFIGURED "${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME}.py")
  configure_file(${TEST_FILE} ${TEST_FILE_CONFIGURED})
  foreach(dependency IN LISTS TEST_DEPENDENCIES)
    configure_file(${dependency} ${CMAKE_CURRENT_BINARY_DIR}/${dependency})
  endforeach(dependency)
  string(REGEX REPLACE "^test_" "${TEST_TYPE}_" TEST_NAME ${TEST_NAME})
  add_test(${TEST_NAME} ${CMAKE_BINARY_DIR}/pypresso ${TEST_FILE_CONFIGURED})
  set_tests_properties(${TEST_NAME} PROPERTIES FIXTURES_REQUIRED
                                               IMPORTLIB_WRAPPER)
  set_tests_properties(${TEST_NAME} PROPERTIES LABELS "${TEST_LABELS}")
  if("gpu" IN_LIST TEST_LABELS)
    set_tests_properties(${TEST_NAME} PROPERTIES RESOURCE_LOCK GPU)
  endif()
endmacro(PYTHON_SCRIPTS_TEST)

add_subdirectory(benchmarks)
add_subdirectory(samples)
add_subdirectory(tutorials)
