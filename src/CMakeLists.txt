# Copyright (C) 2009-2019 The ESPResSo project
# Copyright (C) 2009,2010
#   Max-Planck-Institute for Polymer Research, Theory Group
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

# Target for the unit tests
add_custom_target(
  check_unit_tests COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${TEST_TIMEOUT}
                           $(ARGS) --output-on-failure)

if(WITH_TESTS)
  # Run unit tests on check
  add_dependencies(check check_unit_tests)
endif(WITH_TESTS)

add_subdirectory(utils)
add_subdirectory(shapes)
add_subdirectory(core)
add_subdirectory(config)

add_subdirectory(profiler)

if(SCAFACOS)
  add_subdirectory(scafacos)
endif(SCAFACOS)

if(WITH_SCRIPT_INTERFACE)
  add_subdirectory(script_interface)
endif()

if(WITH_PYTHON)
  add_subdirectory(python)
endif()
