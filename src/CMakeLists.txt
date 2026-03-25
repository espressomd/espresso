#
# Copyright (C) 2009-2026 The ESPResSo project
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

if(ESPRESSO_BUILD_TESTS)
  # Target for the unit tests
  add_custom_target(unit_tests_executables)
  add_custom_target(
    check_unit_tests
    COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${ESPRESSO_TEST_TIMEOUT}
            ${ESPRESSO_CTEST_ARGS} --output-on-failure)

  # Run unit tests on check
  add_dependencies(check check_unit_tests)
  add_dependencies(check_unit_tests unit_tests_executables)
  if(TARGET espresso_packaging_dependencies)
    add_dependencies(espresso_packaging_dependencies unit_tests_executables)
  endif()
endif()

add_subdirectory(utils)
add_subdirectory(shapes)
add_subdirectory(core)
add_subdirectory(config)
add_subdirectory(instrumentation)

if(ESPRESSO_BUILD_WITH_SCAFACOS)
  add_subdirectory(scafacos)
endif()

if(ESPRESSO_BUILD_WITH_WALBERLA)
  add_subdirectory(walberla_bridge)
endif()

if(ESPRESSO_BUILD_WITH_PYTHON)
  add_subdirectory(script_interface)
  add_subdirectory(python)
endif()

add_subdirectory(particle_observables)
