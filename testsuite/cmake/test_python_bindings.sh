#!/usr/bin/env bash

# load bash unit testing library
source BashUnitTests.sh

# test if Python module can be imported
function test_Python() {
  assert_return_code "@PYTHON_EXECUTABLE@" -c "import sys;sys.path.insert(0, '@DESTDIR@/usr/local/lib/python@PYTHON_VERSION_MAJOR@.@PYTHON_VERSION_MINOR@/site-packages');import espressomd"
}

# run tests
run_test_suite

