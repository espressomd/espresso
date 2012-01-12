# Copyright (C) 2009,2010 The ESPResSo project
# Copyright (C) 2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany

# This file is part of ESPResSo.

# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

########################################################################
#Process MPI settings
########################################################################
if(WITH_MPI)
  set(MPI_SOURCE)
else(WITH_MPI)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/mpifake)  
  set(MPI "fake")
  set(MPI_SOURCE mpifake/mpi.c)
endif(WITH_MPI)

########################################################################
# Basic system tests (standard libraries, headers, functions, types)   #
########################################################################

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/acconfig.h.cmakein ${CMAKE_CURRENT_SOURCE_DIR}/acconfig.h)
file(GLOB ESPRESSO_SOURCES *.c)
file(GLOB NOT_ESPRESSO_SOURCES blockfile_test.c scriptsdir.c)
list(REMOVE_ITEM ESPRESSO_SOURCES ${NOT_ESPRESSO_SOURCES})
add_custom_command(OUTPUT scriptsdir-install.c 
  COMMAND ${CMAKE_COMMAND} -E copy_if_different 
  ${CMAKE_CURRENT_SOURCE_DIR}/scriptsdir.c ${CMAKE_CURRENT_BINARY_DIR}/scriptsdir-install.c
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/scriptsdir.c)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES scriptsdir-install.c)
include_directories(${CMAKE_BINARY_DIR}) #for myconfig.h
include_directories(${CMAKE_CURRENT_SOURCE_DIR}) #for all header
set_source_files_properties(scriptsdir-install.c PROPERTIES COMPILE_DEFINITIONS
  ESPRESSO_SCRIPTS_DEFAULT="${CMAKE_PREFIX_PATH}/${SCRIPTDIR}")

add_library(libespresso STATIC ${ESPRESSO_SOURCES} ${CMAKE_BINARY_DIR}/${MYCONFIG} ${MPI_SOURCE})
target_link_libraries(libespresso ${LIBRARIES})
set_target_properties(libespresso PROPERTIES OUTPUT_NAME espresso)

add_executable(Espresso-install scriptsdir-install.c)
target_link_libraries(Espresso-install libespresso)
install(TARGETS Espresso-install RUNTIME DESTINATION bin)
set_target_properties(Espresso-install PROPERTIES OUTPUT_NAME Espresso)
