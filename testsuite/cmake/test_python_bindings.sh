#!/usr/bin/env bash

# load bash unit testing library
source BashUnitTests.sh

# test if Python module can be imported
function test_Python() {
  # test espresso installation via `make install DESTDIR=/some/dir`
  assert_return_code "@CMAKE_BINARY_DIR@/pypresso" -c "import sys;sys.path.insert(0, '@DESTDIR@/@CMAKE_INSTALL_PREFIX@/@Python_SITEARCH@');import espressomd"
  # test espresso installation via `cmake -DCMAKE_INSTALL_PREFIX=/some/dir ..`
  if [ "@CMAKE_INSTALL_PREFIX@" = "/tmp/espresso-unit-tests" ]
  then
    assert_return_code "@CMAKE_BINARY_DIR@/pypresso" -c "import sys;sys.path.insert(0, '@CMAKE_INSTALL_PREFIX@/@Python_SITEARCH@');import espressomd"
  fi
}

# run tests
run_test_suite

