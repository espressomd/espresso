# Detect if compiler version is sufficient for supporting C++11.
# If it is, CMAKE_CXX_FLAGS are modified appropriately and HAVE_CXX11
# is set to a true value.  Else, CMake exits with a fatal error message.
# This currently only works for GCC and Clang compilers.
# In Cmake 3.1+, CMAKE_CXX_STANDARD_REQUIRED should be able to replace
# all the logic below.

if ( DEFINED HAVE_CXX11 )
  return()
endif ()

set(required_gcc_version 4.8)
set(required_clang_version 3.3)
set(required_icc_version 11.0)

# CMAKE_CXX_COMPILER_VERSION may not always be available (e.g. particularly
# for CMakes older than 2.8.10, but use it if it exists.
if ( DEFINED CMAKE_CXX_COMPILER_VERSION )
  if ( CMAKE_CXX_COMPILER_ID STREQUAL "GNU" )
    if ( CMAKE_CXX_COMPILER_VERSION VERSION_LESS ${required_gcc_version} )
      message(FATAL_ERROR "GCC version must be at least "
	"${required_gcc_version} for C++11 support, detected: "
	"${CMAKE_CXX_COMPILER_VERSION}")
    endif ()
  elseif ( CMAKE_CXX_COMPILER_ID STREQUAL
      "Clang" )
    if ( CMAKE_CXX_COMPILER_VERSION VERSION_LESS
	${required_clang_version} )
      message(FATAL_ERROR "Clang version must be at least "
	"${required_clang_version} for C++11 support, detected: "
	"${CMAKE_CXX_COMPILER_VERSION}")
    endif ()
  endif ()

  if( ${CMAKE_CXX_COMPILER} MATCHES "icpc" OR "${CMAKE_CXX_COMPILER_ARG1}" MATCHES "icpc" )
    if( CMAKE_CXX_COMPILER_VERSION VERSION_LESS ${required_icc_version} )
      message(FATAL_ERROR "icpc version must be at least "
	"${required_icc_version} for C++11 support, detected: "
	"${CMAKE_CXX_COMPILER_VERSION}")     
    endif()
  endif()

  if( ${CMAKE_CXX_COMPILER_ID} MATCHES Cray )
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.4)
      message( FATAL_ERROR "Insufficient Cray Compiler Environment version" )
    endif()
    set(HAVE_CXX11 true)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -hstd=c++11")
    return()
  endif()

  set(HAVE_CXX11 true)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  return()
endif ()
