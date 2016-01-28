# Copyright (C) 2009,2010,2011,2012,2015 The ESPResSo project
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

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

message(STATUS "writing ${CMAKE_BINARY_DIR}/myconfig-sample.hpp")
execute_process(COMMAND ${PYTHON_EXECUTABLE}  
                ${CMAKE_CURRENT_SOURCE_DIR}/gen_sampleconfig.py 
                ${CMAKE_SOURCE_DIR}/src/features.def
                OUTPUT_FILE ${CMAKE_BINARY_DIR}/myconfig-sample.hpp
                )



