#
# Copyright (C) 2020-2022 The ESPResSo project
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

function(UTILS_TEST)
  python_scripts_test(${ARGV} TYPE utils) # forward arguments
  set(utils_tests ${utils_tests} ${TEST_FILE_CONFIGURED} PARENT_SCOPE)
endfunction(UTILS_TEST)

utils_test(FILE test_maintainer_CI_jupyter_warnings.py)
utils_test(FILE test_check_myconfig.py)

add_custom_target(
  check_utils
  COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${ESPRESSO_TEST_TIMEOUT}
          -j${ESPRESSO_TEST_NP} ${ESPRESSO_CTEST_ARGS} --output-on-failure)
