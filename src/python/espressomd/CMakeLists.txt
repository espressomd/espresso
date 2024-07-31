#
# Copyright (C) 2013-2022 The ESPResSo project
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

include(espresso_resource_files)

add_custom_target(espressomd)

file(GLOB cython_SRC *.pyx)
file(GLOB cython_HEADER *.pxd)
file(GLOB python_SRC RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" *.py)

espresso_target_resources(espressomd ${python_SRC})

add_subdirectory(io)
add_subdirectory(detail)
add_subdirectory(plugins)

# Install resource files (Python files, text files, etc.)
get_property(ESPRESSOMD_RESOURCE_FILES TARGET espressomd
             PROPERTY EspressoResourceFiles)
foreach(RESOURCE_ABSPATH ${ESPRESSOMD_RESOURCE_FILES})
  cmake_path(RELATIVE_PATH RESOURCE_ABSPATH BASE_DIRECTORY
             ${CMAKE_CURRENT_BINARY_DIR} OUTPUT_VARIABLE RESOURCE_RELPATH)
  cmake_path(GET RESOURCE_RELPATH PARENT_PATH RESOURCE_RELPARENT)
  install(
    FILES "${RESOURCE_ABSPATH}"
    DESTINATION "${ESPRESSO_INSTALL_PYTHON}/espressomd/${RESOURCE_RELPARENT}")
endforeach()

add_library(espresso_pyx_flags INTERFACE)
add_library(espresso::pyx_flags ALIAS espresso_pyx_flags)
target_compile_options(
  espresso_pyx_flags
  INTERFACE
    $<$<NOT:$<CXX_COMPILER_ID:Intel,IntelLLVM>>:-Wno-pedantic>
    $<$<NOT:$<CXX_COMPILER_ID:Intel,IntelLLVM>>:-Wno-cast-qual>
    $<$<NOT:$<CXX_COMPILER_ID:Intel,IntelLLVM>>:-Wno-deprecated-declarations>
    $<$<CXX_COMPILER_ID:Intel,IntelLLVM>:-diag-disable=1224>
    $<$<CXX_COMPILER_ID:GNU>:-Wno-cpp>
    $<$<CXX_COMPILER_ID:GNU>:-Wno-strict-aliasing>
    $<$<CXX_COMPILER_ID:GNU>:-Wno-maybe-uninitialized>
    $<$<CXX_COMPILER_ID:GNU>:-Wno-volatile>
    $<$<CXX_COMPILER_ID:Clang,AppleClang,IntelLLVM>:-Wno-sometimes-uninitialized>
    $<$<CXX_COMPILER_ID:Clang,AppleClang,IntelLLVM>:-Wno-\#warnings>
    $<$<CXX_COMPILER_ID:Clang,AppleClang,IntelLLVM>:-Wno-deprecated-volatile>
    -Wno-missing-field-initializers
    -Wno-unused-variable)

# Configure, compile and install Cython files
foreach(cython_file ${cython_SRC})
  cmake_path(GET cython_file STEM basename)
  cmake_path(RELATIVE_PATH cython_file BASE_DIRECTORY
             ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_VARIABLE relpath)
  cmake_path(GET relpath PARENT_PATH subfolder)
  cmake_path(APPEND CMAKE_CURRENT_BINARY_DIR ${subfolder} "${basename}.cpp"
             OUTPUT_VARIABLE outputpath)
  add_custom_command(
    OUTPUT ${outputpath}
    COMMAND
      ${CYTHON_EXECUTABLE} -3 --cplus --directive embedsignature=True
      --directive binding=True -I ${CMAKE_CURRENT_SOURCE_DIR} -I
      ${CMAKE_CURRENT_BINARY_DIR} ${cython_file} -o ${outputpath}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
    DEPENDS ${cython_file} ${cython_HEADER})
  set(target "espressomd_${basename}")
  add_library(${target} SHARED ${outputpath})
  if(NOT "${subfolder}" STREQUAL "")
    set_target_properties(${target} PROPERTIES PREFIX "${subfolder}/")
  else()
    set_target_properties(${target} PROPERTIES PREFIX "")
  endif()
  set_target_properties(${target} PROPERTIES OUTPUT_NAME ${basename})
  if(APPLE)
    set_target_properties(
      ${target} PROPERTIES SUFFIX ".so" LINK_FLAGS "-undefined dynamic_lookup")
  endif()
  set_target_properties(${target} PROPERTIES CXX_CLANG_TIDY "")
  target_link_libraries(${target} PRIVATE espresso::config espresso::core
                                          espresso::script_interface)
  target_link_libraries(${target} PRIVATE espresso::cpp_flags)
  target_link_libraries(${target} PRIVATE espresso::pyx_flags)
  target_include_directories(
    ${target} SYSTEM PRIVATE ${Python_INCLUDE_DIRS}
                             ${Python_NumPy_INCLUDE_DIRS})
  add_dependencies(espressomd ${target})
  install(TARGETS "${target}"
          LIBRARY DESTINATION "${ESPRESSO_INSTALL_PYTHON}/espressomd")
endforeach()
