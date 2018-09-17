#!/usr/bin/env bash

# load bash unit testing library
source BashUnitTests.sh

# test Sphinx documentation
function test_Sphinx() {
  local OUTDIR="install"
  for filepath in "@CMAKE_BINARY_DIR@/doc/sphinx/modules.rst" \
                  "@CMAKE_BINARY_DIR@/doc/sphinx/html/modules.html" \
                  "@CMAKE_BINARY_DIR@/doc/sphinx/html/index.html" \
                  ;
  do
    assert_file_exists "${filepath}"
  done
}

# run tests
run_test_suite

