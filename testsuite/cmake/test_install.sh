#!/usr/bin/env bash
#
# Copyright (C) 2018-2024 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

# load bash unit testing library
source BashUnitTests.sh

# test installation and Python bindings
function test_install() {
  # check Python files were installed in espressomd
  local -r filepaths=("@ESPRESSO_INSTALL_BINDIR@/pypresso" \
                      "@ESPRESSO_INSTALL_PYTHON@/espressomd/espresso_core.so" \
                      "@ESPRESSO_INSTALL_PYTHON@/espressomd/_init.so" \
                      "@ESPRESSO_INSTALL_PYTHON@/espressomd/__init__.py" \
                      "@ESPRESSO_INSTALL_PYTHON@/espressomd/io/__init__.py" \
                      "@ESPRESSO_INSTALL_PYTHON@/espressomd/io/writer/__init__.py" \
                      "@ESPRESSO_INSTALL_PYTHON@/object_in_fluid/__init__.py" \
                      "@ESPRESSO_INSTALL_PYTHON@/object_in_fluid/oif_classes.py" \
                     )
  for filepath in "${filepaths[@]}"; do
    assert_file_exists "${filepath}"
  done

  # check no Python file was installed outside espressomd
  paths=$(find "@CMAKE_INSTALL_PREFIX@" -path "@ESPRESSO_INSTALL_PYTHON@/espressomd" -prune -o -path "@ESPRESSO_INSTALL_PYTHON@/object_in_fluid" -prune -o \( -name '*.py' -o -name '*.so' \) -print)
  count=$(echo "${paths}" | wc -l)
  assert_string_equal "${paths}" "" "${count} files were installed in the wrong directories:"$'\n'"${paths}"

  # check the espressomd module can be imported from pypresso
  assert_return_code "@CMAKE_INSTALL_FULL_BINDIR@/pypresso" -c "import espressomd"
}

# run tests
run_test_suite
