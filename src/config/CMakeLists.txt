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

add_subdirectory(include/config)
add_subdirectory(src)

configure_file(${CMAKE_SOURCE_DIR}/cmake/espresso_cmake_config.cmakein
               ${CMAKE_CURRENT_BINARY_DIR}/include/config/cmake_config.hpp)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/include/config/config-features.hpp
         ${CMAKE_CURRENT_BINARY_DIR}/src/config-features.cpp
  COMMAND
    ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/gen_featureconfig.py
    ${CMAKE_CURRENT_SOURCE_DIR}/features.def
    ${CMAKE_CURRENT_BINARY_DIR}/include/config/config-features.hpp
    ${CMAKE_CURRENT_BINARY_DIR}/src/config-features.cpp
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/features.def
          ${CMAKE_CURRENT_SOURCE_DIR}/gen_featureconfig.py)
add_custom_target(
  generate_config_features
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/include/config/config-features.hpp
          ${CMAKE_CURRENT_BINARY_DIR}/src/config-features.cpp)

add_custom_target(
  check_myconfig
  COMMAND
    ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_myconfig.py
    ${CMAKE_CXX_COMPILER} ${CMAKE_CURRENT_SOURCE_DIR}/features.def
    ${CMAKE_SOURCE_DIR}/cmake/espresso_cmake_config.cmakein
    ${CMAKE_CURRENT_BINARY_DIR}/include/config/myconfig-final.hpp
    ${CMAKE_CURRENT_BINARY_DIR}/include/config/cmake_config.hpp
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/features.def myconfig)

execute_process(
  COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/gen_sampleconfig.py
          ${CMAKE_CURRENT_SOURCE_DIR}/features.def
  OUTPUT_FILE ${CMAKE_BINARY_DIR}/myconfig-sample.hpp)

add_library(espresso_config SHARED src/config-features.cpp)
add_library(espresso::config ALIAS espresso_config)
add_dependencies(espresso_config myconfig check_myconfig
                 generate_config_features)
install(TARGETS espresso_config
        LIBRARY DESTINATION ${ESPRESSO_INSTALL_PYTHON}/espressomd)
target_include_directories(
  espresso_config PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
                         ${CMAKE_CURRENT_BINARY_DIR}/include)

find_package(Git)
# Parse repository info from git if available Run this at build time to avoid
# rebuilds
add_custom_target(
  version
  COMMAND
    ${CMAKE_COMMAND} -DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}
    -DPROJECT_VERSION=${PROJECT_VERSION} -DGIT_EXECUTABLE=${GIT_EXECUTABLE}
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -P
    ${PROJECT_SOURCE_DIR}/cmake/espresso_version.cmake)
set_property(
  DIRECTORY APPEND
  PROPERTY ADDITIONAL_MAKE_CLEAN_FILES include/config/version.hpp
           include/config/version.hpp.tmp)
add_dependencies(espresso_config version)
