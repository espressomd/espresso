# - Find HDF5
# Find the Hierarchical Data Format 5 (HDF5) includes and library
#
# This module defines
#  HDF5_FOUND
#  HDF5_INCLUDE_DIR
#  HDF5_LIBRARY
#  HDF5_CPP_LIBRARY
#  HDF5_HL_LIBRARY
#  HDF5_LIBRARIES
#

#=============================================================================
# Copyright 2002-2009 Kitware, Inc.
# Copyright 2008-2011 Peter Colberg
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file COPYING-CMAKE-SCRIPTS for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

find_path(HDF5_INCLUDE_DIR hdf5.h
  HINTS
  $ENV{HDF5_DIR} $ENV{HDF5_ROOT} $ENV{HDF5_HOME}
  PATH_SUFFIXES include
)

if(HDF5_USE_STATIC_LIBS)
  set( _HDF5_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
  if(WIN32)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .lib .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
  endif()
endif()

find_library(HDF5_LIBRARY NAMES hdf5
  HINTS
  $ENV{HDF5_DIR} $ENV{HDF5_ROOT} $ENV{HDF5_HOME}
  PATH_SUFFIXES lib64 lib
)

find_library(HDF5_HL_LIBRARY NAMES hdf5_hl
  HINTS
  $ENV{HDF5_DIR} $ENV{HDF5_ROOT} $ENV{HDF5_HOME}
  PATH_SUFFIXES lib64 lib
)

if(HDF5_USE_STATIC_LIBS)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${_HDF5_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})
endif()

if(HDF5_LIBRARY AND HDF5_HL_LIBRARY)
  set(HDF5_LIBRARIES ${HDF5_LIBRARY} ${HDF5_HL_LIBRARY})
endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(HDF5 DEFAULT_MSG
  HDF5_INCLUDE_DIR
  HDF5_LIBRARY
  HDF5_HL_LIBRARY
)

mark_as_advanced(
  HDF5_INCLUDE_DIR
  HDF5_LIBRARY
  HDF5_HL_LIBRARY
)
