#!/usr/bin/env bash

# load bash unit testing library
source BashUnitTests.sh
source test_install.sh

# test espresso installation via `cmake -DCMAKE_INSTALL_PREFIX=/some/dir ..`
function test_install_CMAKE_INSTALL_PREFIX() {
  helper_test_install_common "@CMAKE_INSTALL_PREFIX@"
}

# run tests
run_test_suite

