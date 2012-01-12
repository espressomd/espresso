# Copyright (C) 2011 Olaf Lenz
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
  set(MYCONFIG ENV{ESPRESSO_MYCONFIG})
# test whether MYCONFIG_NAME is a relative/absolute filename
elseif (EXISTS ${MYCONFIG_NAME})
  set(MYCONFIG ${MYCONFIG_NAME})
else()
  # test whether MYCONFIG_NAME is found in the object or source dir
  find_file(MYCONFIG ${MYCONFIG_NAME}
    PATHS ${BINARY_DIR} ${SOURCE_DIR}
    NO_DEFAULT_PATH)
  # use the default if it is not
  if(NOT MYCONFIG)
    set(MYCONFIG ${SOURCE_DIR}/src/myconfig-default.h)
  endif()
endif()

# copy the file to src/myconfig-final.h
execute_process(
  COMMAND ${CMAKE_COMMAND} 
  -E copy_if_different ${MYCONFIG} ${BINARY_DIR}/src/myconfig-final.h
  )

message(STATUS "Config file: ${MYCONFIG}")
