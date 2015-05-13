# Copyright (C) 2009,2010 Christoph Junghans
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
# - Define macro to check inline keyword
#
#  GMX_TEST_INLINE(VARIABLE)
#
#  VARIABLE will be set to the keyword
#
#  Remember to have a cmakedefine for it too...

MACRO(ESP_TEST_INLINE VARIABLE)
    IF(NOT DEFINED TEST_${VARIABLE})

        MESSAGE(STATUS "Checking for inline keyword")

	FOREACH(KEYWORD "static inline" "inline static" "inline" "static")
            IF(NOT TEST_${VARIABLE})
                MESSAGE(STATUS "Checking for inline keyword - try ${KEYWORD}")
                TRY_COMPILE(TEST_${VARIABLE} "${CMAKE_BINARY_DIR}"    
                            "${CMAKE_SOURCE_DIR}/cmake/TestInline.c"
                            COMPILE_DEFINITIONS -D"INLINEDEF=${KEYWORD}" )
                SET(CHK_INLINE_KEYWORD ${KEYWORD})
            ENDIF(NOT TEST_${VARIABLE})
        ENDFOREACH(KEYWORD)
             
        IF(TEST_${VARIABLE})
            SET(${VARIABLE} ${CHK_INLINE_KEYWORD})
            MESSAGE(STATUS "Checking for inline keyword - using ${CHK_INLINE_KEYWORD}")
        ELSE(TEST_${VARIABLE})
	    SET(${VARIABLE} " ")
            MESSAGE(FATAL_ERROR "Checking for inline keyword - none found")
        ENDIF(TEST_${VARIABLE})

    ENDIF(NOT DEFINED TEST_${VARIABLE})        
ENDMACRO(ESP_TEST_INLINE VARIABLE)




