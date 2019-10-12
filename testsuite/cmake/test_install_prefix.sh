#!/usr/bin/env bash
# Copyright (C) 2018-2019 The ESPResSo project
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

# load bash unit testing library
source BashUnitTests.sh
source test_install.sh

# test espresso installation via `cmake -DCMAKE_INSTALL_PREFIX=/some/dir ..`
function test_install_CMAKE_INSTALL_PREFIX() {
  helper_test_install_common "@CMAKE_INSTALL_PREFIX@"
}

# run tests
run_test_suite

