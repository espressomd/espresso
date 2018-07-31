# Copyright (C) 2009,2010,2011,2012,2015,2016 The ESPResSo project
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
cmake_minimum_required(VERSION 2.8)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
if(SCAFACOS)
  include_directories(${SCAFACOS_INCLUDE_DIRS})
endif(SCAFACOS)

configure_file(${CMAKE_SOURCE_DIR}/cmake/cmake_config.cmakein
               cmake_config.hpp
)

add_custom_target(check_myconfig
  COMMAND ${PYTHON_EXECUTABLE}
  ${CMAKE_CURRENT_SOURCE_DIR}/check_myconfig.py
  ${CMAKE_CXX_COMPILER}
  ${CMAKE_CURRENT_SOURCE_DIR}/features.def
  ${CMAKE_CURRENT_BINARY_DIR}/core/myconfig-final.hpp
  cmake_config.hpp
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/features.def
          myconfig
  )

add_subdirectory(core)

if(WITH_SCRIPT_INTERFACE)
  add_subdirectory(script_interface)
endif()


if(WITH_PYTHON)
  add_subdirectory(python)
endif()

