#
# Copyright (C) 2019-2022 The ESPResSo project
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

function(SCRIPTS_TEST)
  espresso_importlib_based_test(${ARGV} TYPE scripts)
endfunction(SCRIPTS_TEST)

scripts_test(FILE test_notebook_links.py)
scripts_test(FILE test_check_myconfig.py)

add_custom_target(
  check_scripts
  COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${ESPRESSO_TEST_TIMEOUT}
          -j${ESPRESSO_TEST_NP} ${ESPRESSO_CTEST_ARGS} --output-on-failure)
