#!/usr/bin/env bash

# load bash unit testing library
source BashUnitTests.sh

# test binaries and libraries are correclty installed
function test_install() {
  local filepaths=("@DESTDIR@/usr/local/bin/ipypresso" \
                   "@DESTDIR@/usr/local/bin/pypresso" \
                   "@DESTDIR@/usr/local/lib/libEspressoCore.so" \
                   "@DESTDIR@/usr/local/lib/libEspressoCore.so"
                  )
  if [ ! -z "@PYTHON_VERSION_STRING@" ]
  then
    local python_dir="@DESTDIR@/usr/local/lib/python@PYTHON_VERSION_MAJOR@.@PYTHON_VERSION_MINOR@"
    filepaths+=("${python_dir}/site-packages/espressomd/_init.so" \
                "${python_dir}/site-packages/espressomd/__init__.py"
               )
  fi
  for filepath in ${filepaths[@]}
  do
    assert_file_exists "${filepath}"
  done
}

# run tests
run_test_suite

