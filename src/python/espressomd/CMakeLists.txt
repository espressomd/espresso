#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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

file(GLOB cython_SRC *.pyx)
file(GLOB cython_HEADER *.pxd)
file(GLOB cython_AUX *.py)
file(COPY ${cython_AUX} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_custom_command(OUTPUT gen_pxiconfig.cpp
                   COMMAND ${PYTHON_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/gen_pxiconfig.py
                   ${CMAKE_SOURCE_DIR}/src/features.def ${CMAKE_CURRENT_BINARY_DIR}/gen_pxiconfig.cpp
                   DEPENDS ${CMAKE_SOURCE_DIR}/src/features.def
                  )

include_directories(${CMAKE_BINARY_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/src/core)
include_directories(${CMAKE_BINARY_DIR}/src/core)

add_executable(gen_pxiconfig gen_pxiconfig.cpp)
add_dependencies(gen_pxiconfig EspressoConfig)

add_custom_command(OUTPUT _espresso/myconfig.pxi
		   COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/_espresso
                   COMMAND ${CMAKE_CURRENT_BINARY_DIR}/gen_pxiconfig > ${CMAKE_CURRENT_BINARY_DIR}/_espresso/myconfig.pxi
                   DEPENDS gen_pxiconfig
                   )

include_directories(${PYTHON_INCLUDE_DIRS})
include_directories(${NUMPY_INCLUDE_DIR})

add_custom_target(espressomd)

foreach(cython_file ${cython_SRC})
  get_filename_component(basename ${cython_file} NAME_WE)

  if("${basename}" STREQUAL "")
      message( FATAL_ERROR "Internal error empty basename of file ${cython_file}")
  else()
      list(APPEND cython_generated_SRC "${basename}.cpp")
      add_custom_command(OUTPUT ${basename}.cpp
                         COMMAND ${CYTHON_EXECUTABLE}  --cplus
                         -I ${CMAKE_CURRENT_BINARY_DIR}/_espresso
                         -I ${CMAKE_CURRENT_SOURCE_DIR}
                         -I ${CMAKE_CURRENT_SOURCE_DIR}/_espresso
                         ${cython_file}
                         -o ${CMAKE_CURRENT_BINARY_DIR}/${basename}.cpp
                         DEPENDS _espresso/myconfig.pxi ${cython_file}
                         )

      add_library(${basename} SHARED ${basename}.cpp)
      set_target_properties(${basename} PROPERTIES PREFIX "")
      set_target_properties(${basename} PROPERTIES MACOSX_RPATH TRUE)
      if (APPLE)
         set_target_properties(${basename} PROPERTIES SUFFIX ".so")
      endif()
      target_link_libraries(${basename} EspressoCore ${LIBRARIES})
      add_dependencies(espressomd ${basename})
  endif()
endforeach()

