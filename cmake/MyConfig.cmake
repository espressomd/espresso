# Copyright (C) 2011 Olaf Lenz
# Copyright (C) 2015 The ESPResSo Project
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
# This cmake script is used to find out what myconfig header to use.
# It needs the variables

# heed the environment variable "ESPRESSO_MYCONFIG"
if (ENV{ESPRESSO_MYCONFIG})
  set(MYCONFIG_FILE ENV{ESPRESSO_MYCONFIG})
else()
  # test whether MYCONFIG_NAME is found in the object or source dir
  find_file(MYCONFIG_FILE 
    NAMES ${MYCONFIG_NAME}
    PATHS ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}
    NO_DEFAULT_PATH)
  # use the default if it is not
  if(NOT MYCONFIG_FILE)
    set(MYCONFIG_FILE ${CMAKE_SOURCE_DIR}/src/core/myconfig-default.hpp)
  endif()
endif()
add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/src/core/myconfig-final.hpp
  COMMAND ${CMAKE_COMMAND} -E copy ${MYCONFIG_FILE} ${CMAKE_BINARY_DIR}/src/core/myconfig-final.hpp
  DEPENDS ${MYCONFIG_FILE})
add_custom_target(myconfig DEPENDS ${CMAKE_BINARY_DIR}/src/core/myconfig-final.hpp)
message(STATUS "Config file: ${MYCONFIG_FILE}")
# Clear variable, otherwise cmake must be run by hand to detect myconfig
# Also prevents find_file from skipping when variable is already set
unset(MYCONFIG_FILE CACHE)
