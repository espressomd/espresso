# Copyright (C) 2009,2010,2011,2012 The ESPResSo project
# Copyright (C) 2009,2010 
#  Max-Planck-Institute for Polymer Research, Theory Group
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

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING 
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

########################################################################
# User input options                                                   #
########################################################################
option(WITH_MPI    "Build a parallel (message-passing) version of ESPResSo" ON)
option(WITH_TK     "Build with tk support" OFF)
set(MYCONFIG "myconfig.h" CACHE STRING "default name of the local config file")
if (NOT DEFINED DATA)
  set(DATA "share/espresso")
endif(NOT DEFINED DATA)
set(SCRIPTDIR ${DATA}/scripts)

# TODO: not portable
enable_language(C)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99") 

find_package(FFTW3)
if (FFTW3_FOUND)
    include_directories(${FFTW3_INCLUDE_DIR})
    list(APPEND LIBRARIES ${FFTW3_LIBRARIES})
    set(FFTW 3)
else(FFTW3_FOUND)
    message(STATUS "FFTW3 not found, disabling it.")
endif(FFTW3_FOUND)

find_package(TCL REQUIRED)
include_directories(${TCL_INCLUDE_PATH})
list(APPEND LIBRARIES ${TCL_LIBRARY})

########################################################################
# Fix stupid flags on MSVC
########################################################################
IF(CMAKE_GENERATOR MATCHES "Visual Studio")
    STRING(REPLACE /MD /MT CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})
    SET(CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE} CACHE STRING "" FORCE)
    STRING(REPLACE /MD /MT CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG})
    SET(CMAKE_C_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG} CACHE STRING "" FORCE)   
ENDIF(CMAKE_GENERATOR MATCHES "Visual Studio")

########################################################################
# Process MPI settings
########################################################################
if(WITH_MPI)
  find_package(MPI)
  if(MPI_FOUND)
    list(APPEND CMAKE_C_FLAGS ${MPI_COMPILE_FLAGS})
    list(APPEND CMAKE_EXE_LINKER_FLAGS ${MPI_LINK_FLAGS})
    include_directories(${MPI_INCLUDE_PATH})
    list(APPEND LIBRARIES ${MPI_LIBRARIES})
    set(MPI "external")
  else(MPI_FOUND)
    message(FATAL_ERROR "MPI support requested, but no MPI compiler found.")
  endif(MPI_FOUND)
endif(WITH_MPI)

include(EspTestInline)
esp_test_inline(MDINLINE)
message(STATUS "Inline keyword: ${MDINLINE}")

########################################################################
# Basic system tests (standard libraries, headers, functions, types)   #
########################################################################
message(STATUS "Using config file: ${MYCONFIG}")
if (NOT EXISTS ${CMAKE_BINARY_DIR}/${MYCONFIG})
  message(FATAL_ERROR "Could not find ${MYCONFIG} in ${CMAKE_BINARY_DIR}")
endif (NOT EXISTS ${CMAKE_BINARY_DIR}/${MYCONFIG})

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

add_subdirectory(src)
add_subdirectory(scripts)
add_executable(Espresso src/scriptsdir.c)
set_source_files_properties(src/scriptsdir.c PROPERTIES COMPILE_DEFINITIONS
  ESPRESSO_SCRIPTS_DEFAULT="${CMAKE_SOURCE_DIR}/scripts")
target_link_libraries(Espresso libespresso)
enable_testing()
add_subdirectory(testsuite)

