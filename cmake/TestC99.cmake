# Copyright (C) 2012 Olaf Lenz, Florian Fahrenberger
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
MACRO(TEST_C99 OUTPUT_VARIABLE)
  ENABLE_LANGUAGE(C)
  SET(ORIG_CMAKE_C_FLAGS ${CMAKE_C_FLAGS})

  FOREACH(TEST_C99_FLAG '' -std=gnu99 -std=c99 -c99 -AC99 -xc99=all -qlanglvl=extc99)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${TEST_C99_FLAG}")
    TRY_COMPILE(FOUND_C99_FLAG ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY} ${CMAKE_MODULE_PATH}/TestC99.c)
    SET(CMAKE_C_FLAGS "${ORIG_CMAKE_C_FLAGS}")
    IF(FOUND_C99_FLAG)
      SET(${OUTPUT_VARIABLE} ${TEST_C99_FLAG})
      BREAK()
    ENDIF(FOUND_C99_FLAG)
  ENDFOREACH(TEST_C99_FLAG)
  IF(NOT FOUND_C99_FLAG)
    SET(${OUTPUT_VARIABLE} C99-NOTFOUND)
  ENDIF(NOT FOUND_C99_FLAG)
ENDMACRO(TEST_C99)
