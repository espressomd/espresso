#!/usr/bin/env bash

# test espresso installation
function helper_test_install_common() {
  local root=$1
  local filepaths=("${root}/bin/pypresso" \
                   "${root}/@Python_SITEARCH@/espressomd/EspressoCore.so" \
                  )
  if [ "@TESTING_PYTHON@" = "TRUE" ]
  then
    filepaths+=("${root}/@Python_SITEARCH@/espressomd/_init.so" \
                "${root}/@Python_SITEARCH@/espressomd/__init__.py"
               )
  fi
  for filepath in ${filepaths[@]}
  do
    assert_file_exists "${filepath}"
  done
}

