#
# Copyright (C) 2015-2022 The ESPResSo project
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

# Configure pypresso for build dir
set(ESPRESSO_PYTHON_DIR ${CMAKE_CURRENT_BINARY_DIR})

set(ESPRESSO_PYTHON_FRONTEND ${Python_EXECUTABLE})
configure_file(pypresso.cmakein ${CMAKE_BINARY_DIR}/pypresso @ONLY)

if(IPYTHON_EXECUTABLE)
  set(ESPRESSO_PYTHON_FRONTEND ${IPYTHON_EXECUTABLE})
  configure_file(pypresso.cmakein ${CMAKE_BINARY_DIR}/ipypresso @ONLY)
endif()

# Configure pypresso for install dir
set(ESPRESSO_PYTHON_DIR ${ESPRESSO_INSTALL_PYTHON})

set(ESPRESSO_PYTHON_FRONTEND ${Python_EXECUTABLE})
configure_file(pypresso.cmakein ${CMAKE_CURRENT_BINARY_DIR}/pypresso @ONLY)

if(IPYTHON_EXECUTABLE)
  set(ESPRESSO_PYTHON_FRONTEND ${IPYTHON_EXECUTABLE})
  configure_file(pypresso.cmakein ${CMAKE_CURRENT_BINARY_DIR}/ipypresso @ONLY)
endif()

option(INSTALL_PYPRESSO
       "Install pypresso script, not needed when ESPResSo is installed in /usr"
       ON)
if(INSTALL_PYPRESSO)
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/pypresso
          DESTINATION ${CMAKE_INSTALL_BINDIR})

  if(IPYTHON_EXECUTABLE)
    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ipypresso
            DESTINATION ${CMAKE_INSTALL_BINDIR})
  endif()
endif(INSTALL_PYPRESSO)

# Make a target for pypresso, so that things can depend on it
add_custom_target(pypresso ALL)
add_dependencies(pypresso espressomd object_in_fluid)

add_subdirectory(espressomd)
add_subdirectory(object_in_fluid)
