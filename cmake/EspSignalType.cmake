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
# - Define macro to check return type of signals (int/void)
#
#  ESP_TEST_RETSIGTYPE(VARIABLE)
#
#  VARIABLE will be set to the return type of signals - "int" or "void"
#
#  Remember to have a cmakedefine for it too...
#
#


MACRO(ESP_TEST_RETSIGTYPE VARIABLE)
    IF(NOT DEFINED ${VARIABLE})

        MESSAGE(STATUS "Checking for return type of signals")

	# First check without any special flags
        TRY_COMPILE(RETSIGTYPE_INT_OK "${CMAKE_BINARY_DIR}"    
                    "${CMAKE_SOURCE_DIR}/cmake/TestRetSigType.c")

        if(RETSIGTYPE_INT_OK)
	    MESSAGE(STATUS "Checking for return type of signals - int")			
            set(${VARIABLE} "int" CACHE INTERNAL "Result of test for signal return type" FORCE)
        else(RETSIGTYPE_INT_OK)
            MESSAGE(STATUS "Checking for return type of signals - void")
      	    set(${VARIABLE} "void" CACHE INTERNAL "Result of test for signal return type" FORCE)
      	endif(RETSIGTYPE_INT_OK)
        
    ENDIF(NOT DEFINED ${VARIABLE})
ENDMACRO(ESP_TEST_RETSIGTYPE VARIABLE)



