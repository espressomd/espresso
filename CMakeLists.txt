# Copyright (C) 2009,2010,2011,2012,2015,2016 The ESPResSo project
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

cmake_minimum_required(VERSION 3.0)

project(ESPResSo)

enable_language(CXX)

set(PROJECT_VERSION "4.0-dev")
string(REGEX REPLACE "^([1-9]+)\\.[0-9]+.*$" "\\1" SOVERSION "${PROJECT_VERSION}")
if (NOT ${SOVERSION} MATCHES "^[1-9]+$")
  message(FATAL_ERROR "Could not determind SOVERSION from ${PROJECT_VERSION}")
endif (NOT ${SOVERSION} MATCHES "^[1-9]+$")

######################################################################
# CMake internal vars
######################################################################

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)


# On Mac OS X, first look for other packages, then frameworks
set(CMAKE_FIND_FRAMEWORK LAST)

######################################################################
# User input options
######################################################################

option(WITH_PYTHON "Build python interface"  ON)
option(WITH_TCL    "Build tcl interface"     ON)
option(WITH_CUDA   "Build with GPU support"  ON)
option(WITH_HDF5   "Build with HDF5 support" ON)
option(WITH_TESTS  "Enable tests"            ON)
option(WITH_SCAFACOS "Build with Scafacos support" ON)
option(WITH_VALGRIND_INSTRUMENTATION "Build with valgrind instrumentation markers" OFF)

# choose the name of the config file
set(MYCONFIG_NAME "myconfig.hpp"
  CACHE STRING "Default name of the local config file")

# Check which config file to use
include(MyConfig)

######################################################################
# Pretty function
######################################################################

include(CheckCXXSourceCompiles)

foreach(func_name __PRETTY_FUNCTION__ __FUNCTION__)
   check_cxx_source_compiles("
     #include <string>
     int main() { std::string(${func_name}); }
     " result${func_name})
   if(result${func_name})
     set(__PRETTYFUNC__ ${func_name})
     break()
   endif(result${func_name})
endforeach(func_name __PRETTY_FUNCTION__ __FUNCTION__)

######################################################################
# Libraries
######################################################################

if (WITH_CUDA)
  find_package(CUDA QUIET)
  if(CUDA_FOUND)
    list(APPEND LIBRARIES ${CUDA_CUFFT_LIBRARIES})
    list(APPEND LIBRARIES ${CUDA_LIBRARIES})
    set(CUDA 1)

    if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
      set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -g -G")
    else()
      set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -O3")
    endif()
    set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -gencode=arch=compute_20,code=sm_20 -gencode=arch=compute_30,code=sm_30 -gencode=arch=compute_50,code=sm_50")
    list(APPEND CUDA_NVCC_FLAGS "-std=c++11")
    SET(CUDA_PROPAGATE_HOST_FLAGS OFF)

    if (APPLE)
      set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -Xcompiler -isysroot -Xcompiler ${CMAKE_OSX_SYSROOT}")
    endif()
  endif(CUDA_FOUND)
endif(WITH_CUDA)

find_package(PythonInterp)

if (WITH_PYTHON)
  find_package(Cython 0.23 REQUIRED)
  execute_process(COMMAND ${PYTHON_EXECUTABLE}
                -c "import distutils.sysconfig as cg; print(cg.get_python_inc())"
                OUTPUT_VARIABLE PYTHON_INCLUDE_DIRS OUTPUT_STRIP_TRAILING_WHITESPACE)
  find_package(NumPy REQUIRED)
  execute_process(COMMAND ${PYTHON_EXECUTABLE} 
                -c "import distutils.sysconfig as cg; print(cg.get_python_lib(1,0,prefix='${CMAKE_INSTALL_EXEC_PREFIX}'))"
                OUTPUT_VARIABLE PYTHON_INSTDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(PYTHON_VERSION_MAJOR EQUAL 2)
    set(CYTHON_FLAGS "-2" CACHE STRING "Flags used by the Cython compiler during all build types.")
  else()
    set(CYTHON_FLAGS "-3" CACHE STRING "Flags used by the Cython compiler during all build types.")
  endif()
endif(WITH_PYTHON)

find_package(FFTW3)
if (FFTW3_FOUND)
  include_directories(${FFTW3_INCLUDE_DIR})
  list(APPEND LIBRARIES ${FFTW3_LIBRARIES})
  set(FFTW 3)
endif(FFTW3_FOUND)

if (WITH_TCL)
  find_package(TCL)

  if(NOT TCL_FOUND)
	  message(FATAL_ERROR "TCL not found.")
  endif()

  include_directories(${TCL_INCLUDE_PATH})
  list(APPEND LIBRARIES ${TCL_LIBRARY})
endif()

if (WITH_HDF5)
find_package(HDF5)
if(HDF5_FOUND)
  set(H5MD 1)
  list(APPEND LIBRARIES ${HDF5_LIBRARIES})
  include_directories(${HDF5_INCLUDE_DIRS})
endif(HDF5_FOUND)
endif(WITH_HDF5)

if(WITH_SCAFACOS)
  find_package(PkgConfig)
  pkg_check_modules(SCAFACOS scafacos)
  if(SCAFACOS_FOUND)
    set(SCAFACOS 1)
  endif(SCAFACOS_FOUND)
endif(WITH_SCAFACOS)

if(WITH_VALGRIND_INSTRUMENTATION)
  find_package(PkgConfig)
  pkg_check_modules(VALGRIND valgrind)
  if(VALGRIND_FOUND)
    set(VALGRIND_INSTRUMENTATION 1)
    message(STATUS ${VALGRIND_INCLUDE_DIRS})
    include_directories(${VALGRIND_INCLUDE_DIRS})
  endif(VALGRIND_FOUND)
endif(WITH_VALGRIND_INSTRUMENTATION)

include(RequireCXX11)

#######################################################################
# Boost
#######################################################################

find_package(Boost REQUIRED mpi serialization)
include_directories(${Boost_INCLUDE_DIRS})
list(APPEND LIBRARIES ${Boost_LIBRARIES})

#######################################################################
# Testing 
#######################################################################

if(WITH_TESTS)
  enable_testing()
  find_package(Boost COMPONENTS unit_test_framework) 
  if(Boost_UNIT_TEST_FRAMEWORK_FOUND)
    set(WITH_UNIT_TESTS ON)
  endif(Boost_UNIT_TEST_FRAMEWORK_FOUND)
  add_custom_target(check)
  add_subdirectory(testsuite)
  list(APPEND LIBRARIES ${Boost_LIBRARIES})
endif(WITH_TESTS)

#######################################################################
# Process MPI settings
#######################################################################

find_package(MPI REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MPI_LINK_FLAGS}")
include_directories(${MPI_INCLUDE_PATH})
list(APPEND LIBRARIES ${MPI_LIBRARIES})

#######################################################################
# Paths
#######################################################################

if (NOT DEFINED DATA)
  set(DATA "share/espresso")
endif(NOT DEFINED DATA)

if (NOT DEFINED LIBDIR)
  set(LIBDIR "lib")
endif(NOT DEFINED LIBDIR)

if (NOT DEFINED BINDIR)
  set(BINDIR "bin")
endif(NOT DEFINED BINDIR)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${LIBDIR}")

set(SCRIPTDIR ${DATA}/scripts)

#######################################################################
# Flags
#######################################################################

if(APPLE)
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -undefined dynamic_lookup")
endif()

#######################################################################
# Subdirectories
#######################################################################

add_subdirectory(scripts)
add_subdirectory(doc)
add_subdirectory(src)
add_subdirectory(config)

#######################################################################
# Feature summary
#######################################################################

include(FeatureSummary)
feature_summary(WHAT ALL)
