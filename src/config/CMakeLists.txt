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

configure_file(${CMAKE_SOURCE_DIR}/cmake/cmake_config.cmakein cmake_config.hpp)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/config-features.hpp
         ${CMAKE_CURRENT_BINARY_DIR}/config-features.cpp
  COMMAND
    ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/gen_featureconfig.py
    ${CMAKE_CURRENT_SOURCE_DIR}/features.def
    ${CMAKE_CURRENT_BINARY_DIR}/config-features.hpp
    ${CMAKE_CURRENT_BINARY_DIR}/config-features.cpp
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/features.def
          ${CMAKE_CURRENT_SOURCE_DIR}/gen_featureconfig.py)
add_custom_target(
  generate_config_features
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/config-features.hpp
          ${CMAKE_CURRENT_BINARY_DIR}/config-features.cpp)

add_custom_target(
  check_myconfig
  COMMAND
    ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_myconfig.py
    ${CMAKE_CXX_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/features.def
    ${CMAKE_CURRENT_BINARY_DIR}/myconfig-final.hpp cmake_config.hpp
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/features.def myconfig)

execute_process(
  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/gen_sampleconfig.py
          ${CMAKE_CURRENT_SOURCE_DIR}/features.def
  OUTPUT_FILE ${CMAKE_BINARY_DIR}/myconfig-sample.hpp)

add_library(Espresso_config SHARED config-features.cpp)
add_library(Espresso::config ALIAS Espresso_config)
add_dependencies(Espresso_config myconfig check_myconfig
                 generate_config_features)
install(TARGETS Espresso_config
        LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)
target_include_directories(Espresso_config PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
                                                  ${CMAKE_CURRENT_BINARY_DIR})

find_package(Git)
# Parse repository info from git if available Run this at build time to avoid
# rebuilds
add_custom_target(
  version
  COMMAND
    ${CMAKE_COMMAND} -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
    -DPROJECT_VERSION=${PROJECT_VERSION} -DGIT_EXECUTABLE=${GIT_EXECUTABLE}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -P
    ${PROJECT_SOURCE_DIR}/cmake/version.cmake)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES version.hpp
                                       version.hpp.tmp)
add_dependencies(Espresso_config version)
