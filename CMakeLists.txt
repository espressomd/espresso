# Copyright (C) 2009,2010,2011,2012 The ESPResSo project
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
cmake_minimum_required(VERSION 2.6)

project(ESPResSo)
set(PROJECT_VERSION "3.0.2")

enable_language(C)

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

######################################################################
# User input options                                                   
######################################################################
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

option(WITH_MPI    "Build a parallel (message-passing) version of ESPResSo" ON)
option(WITH_TK     "Build with tk support" OFF)

# choose the name of the config file
set(ESPRESSO_MYCONFIG_NAME "myconfig.h"
  CACHE STRING "Default name of the local config file")

###################################################################### 
# Compiler characteristics
######################################################################
include(TestC99)
test_c99(C99_FLAG)
if(C99_FLAG)
  message(STATUS "Compiler option for C99: ${C99_FLAG}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${C99_FLAG}")
else(C99_FLAG)
  message(FATAL_ERROR "Compiling ESPResSo requires a C99 compiler.")
endif(C99_FLAG)

include(EspTestInline)
esp_test_inline(MDINLINE)
message(STATUS "Inline keyword: ${MDINLINE}")

######################################################################
# Fix stupid flags on MSVC
######################################################################
IF(CMAKE_GENERATOR MATCHES "Visual Studio")
    STRING(REPLACE /MD /MT CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})
    SET(CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE} CACHE STRING "" FORCE)
    STRING(REPLACE /MD /MT CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG})
    SET(CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG} CACHE STRING "" FORCE)   
ENDIF(CMAKE_GENERATOR MATCHES "Visual Studio")

######################################################################
# Libraries
######################################################################
#
find_package(FFTW3)
if (FFTW3_FOUND)
    include_directories(${FFTW3_INCLUDE_DIR})
    list(APPEND LIBRARIES ${FFTW3_LIBRARIES})
    set(FFTW 3)
else(FFTW3_FOUND)
    message(STATUS "FFTW3 not found, diabling it.")
endif(FFTW3_FOUND)

find_package(TCL REQUIRED)
include_directories(${TCL_INCLUDE_PATH})
list(APPEND LIBRARIES ${TCL_LIBRARY})

#######################################################################
# Process MPI settings
#######################################################################
if(WITH_MPI)
  find_package(MPI)
  if(MPI_FOUND)
    list(APPEND CMAKE_C_FLAGS ${MPI_COMPILE_FLAGS})
    list(APPEND CMAKE_EXE_LINKER_FLAGS ${MPI_LINK_FLAGS})
    include_directories(${MPI_INCLUDE_PATH})
    list(APPEND LIBRARIES ${MPI_LIBRARIES})
    set(MPI_SOURCE)
  else(MPI_FOUND)
    message(FATAL_ERROR "MPI support requested, but no MPI compiler found.")
  endif(MPI_FOUND)
else(WITH_MPI)
  include_directories(${CMAKE_SOURCE_DIR}/src/mpifake)  
  set(MPI_SOURCE src/mpifake/mpi.c)
endif(WITH_MPI)

#######################################################################
# Basic system tests (standard libraries, headers, functions, types)
#######################################################################

include(CheckIncludeFiles)
check_include_files("arpa/inet.h;ctype.h;errno.h;fcntl.h;float.h;inttypes.h;limits.h;math.h;netdb.h;netinet/in.h;signal.h;stddef.h;stdio.h;stdlib.h;string.h;strings.h;sys/file.h;sys/resource.h;sys/socket.h;sys/time.h;sys/types.h;time.h;unistd.h" HAVE_HEADERS)
if (NOT HAVE_HEADERS)
  message(FATAL_ERROR "There is a least one header missing.")
endif(NOT HAVE_HEADERS)

include(CheckFunctionExists)
check_function_exists("atexit;floor;getcwd;gethostbyname;memset;pow;rint;select;socket;sqrt"       HAVE_FUNCS)

include(CheckLibraryExists)
check_library_exists(m sqrt "" HAVE_LIBM)

include(CheckTypeSize)
check_type_size(size_t size_t)
if(NOT size_t)
    set(size_t int)
else(NOT size_t)
    set(size_t 0)
endif(NOT size_t)

include(EspSignalType)
ESP_TEST_RETSIGTYPE(RETSIGTYPE)

set(MATH_LIBRARIES "m" CACHE STRING "math library")
mark_as_advanced( MATH_LIBRARIES )
include(CheckLibraryExists)
foreach(FUNC sqrt)
  check_library_exists(${MATH_LIBRARIES} ${FUNC} "" FOUND_${FUNC}_${MATH_LIBRARIES})
  if(NOT FOUND_${FUNC}_${MATH_LIBRARIES})
    message(FATAL_ERROR "Could not find needed math function - ${FUNC}")
  endif(NOT FOUND_${FUNC}_${MATH_LIBRARIES})
endforeach(FUNC)
list(APPEND LIBRARIES ${MATH_LIBRARIES})


######################################################################
# Specify install locations and which subdirectories to process
######################################################################
if (NOT DEFINED DATA)
  set(DATA "share/espresso")
endif(NOT DEFINED DATA)
set(SCRIPTDIR ${DATA}/scripts)

add_subdirectory(scripts)

######################################################################
# Compile binary
######################################################################
# generate acconfig.h
configure_file(${CMAKE_SOURCE_DIR}/src/acconfig.h.in ${CMAKE_BINARY_DIR}/src/acconfig.h)

# include dirs
include_directories(${CMAKE_SOURCE_DIR}/src ${CMAKE_BINARY_DIR}/src)

# collect source files
file(GLOB ESPRESSO_SOURCES src/*.c)
file(GLOB NOT_ESPRESSO_SOURCES src/*_test.c)
file(GLOB ESPRESSO_TCL_SOURCES src/tcl/*.c)
list(REMOVE_ITEM ESPRESSO_SOURCES ${NOT_ESPRESSO_SOURCES})

# compile executable
set_source_files_properties(src/scriptsdir.c PROPERTIES 
  COMPILE_DEFINITIONS
  ESPRESSO_SCRIPTS_DEFAULT="${CMAKE_BINARY_DIR}/scripts")
add_executable(Espresso 
  ${ESPRESSO_SOURCES} 
  ${ESPRESSO_TCL_SOURCES} 
  ${MPI_SOURCE}
  )
target_link_libraries(Espresso ${LIBRARIES})
add_dependencies(Espresso scripts)

# handle myconfig
add_dependencies(Espresso myconfig-header scripts)

install(TARGETS Espresso RUNTIME DESTINATION bin)

enable_testing()
add_subdirectory(testsuite)

######################################################################
# Handling of myconfig
###################################################################### 
# Execute FindMyConfig.cmake whenever ESPResSo is built
add_custom_target(myconfig-header
  COMMAND ${CMAKE_COMMAND} 
  -D SOURCE_DIR=${CMAKE_SOURCE_DIR}
  -D BINARY_DIR=${CMAKE_BINARY_DIR}
  -D MYCONFIG_NAME=${ESPRESSO_MYCONFIG_NAME}
  -P ${CMAKE_SOURCE_DIR}/cmake/FindMyConfig.cmake
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMENT "Determining config file."
  )

