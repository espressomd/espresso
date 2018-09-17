#!/usr/bin/env bash

# load bash unit testing library
source BashUnitTests.sh

# test Doxygen documentation
function test_Doxygen() {
  assert_non_zero "$(ls @CMAKE_BINARY_DIR@/doc/doxygen/html/*_8hpp_source.html    | wc -l)" "No output for .hpp sources"
  assert_non_zero "$(ls @CMAKE_BINARY_DIR@/doc/doxygen/html/*_8hpp.html           | wc -l)" "No output for C++ headers"
  assert_non_zero "$(ls @CMAKE_BINARY_DIR@/doc/doxygen/html/*_8cpp_source.html    | wc -l)" "No output for .cpp sources"
  assert_non_zero "$(ls @CMAKE_BINARY_DIR@/doc/doxygen/html/*_8cpp.html           | wc -l)" "No output for C++ code"
  assert_non_zero "$(ls @CMAKE_BINARY_DIR@/doc/doxygen/html/*_8cu_source.html     | wc -l)" "No output for .cu sources"
  assert_non_zero "$(ls @CMAKE_BINARY_DIR@/doc/doxygen/html/*_8cu.html            | wc -l)" "No output for C++ CUDA"
  assert_non_zero "$(ls @CMAKE_BINARY_DIR@/doc/doxygen/html/struct*.html          | wc -l)" "No output for C++ structs"
  assert_non_zero "$(ls @CMAKE_BINARY_DIR@/doc/doxygen/html/class*.html           | wc -l)" "No output for C++ classes"
  assert_non_zero "$(ls @CMAKE_BINARY_DIR@/doc/doxygen/html/namespace*.html       | wc -l)" "No output for C++ namespaces"
  assert_non_zero "$(ls @CMAKE_BINARY_DIR@/doc/doxygen/html/functions_func_*.html | wc -l)" "No output for C++ functions/methods"
  assert_non_zero "$(ls @CMAKE_BINARY_DIR@/doc/doxygen/html/functions_vars_*.html | wc -l)" "No output for C++ variables/members"
}

# run tests
run_test_suite

