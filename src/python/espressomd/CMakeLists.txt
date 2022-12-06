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

add_custom_command(
  OUTPUT gen_pxiconfig.cpp
  COMMAND
    ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/gen_pxiconfig.py
    ${CMAKE_SOURCE_DIR}/src/config/features.def
    ${CMAKE_CURRENT_BINARY_DIR}/gen_pxiconfig.cpp
  DEPENDS ${CMAKE_SOURCE_DIR}/src/config/features.def)

add_executable(gen_pxiconfig gen_pxiconfig.cpp)
target_link_libraries(gen_pxiconfig espresso::config)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/myconfig.pxi
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/gen_pxiconfig >
          ${CMAKE_CURRENT_BINARY_DIR}/myconfig.pxi DEPENDS gen_pxiconfig)

add_custom_target(espressomd)

# Make the cython_SRC, cython_HEADER and cython_AUX a cached variable to be able
# to extend it in the subdirectories.
file(GLOB cython_SRC *.pyx)
set(cython_SRC "${cython_SRC}" CACHE INTERNAL "cython_SRC")
file(GLOB cython_HEADER *.pxd)
set(cython_HEADER "${cython_HEADER}" CACHE INTERNAL "cython_HEADER")
file(GLOB cython_AUX *.py)
set(cython_AUX "${cython_AUX}" CACHE INTERNAL "cython_AUX")

add_subdirectory(io)

list(REMOVE_DUPLICATES cython_SRC)

add_library(espresso_pyx_flags INTERFACE)
add_library(espresso::pyx_flags ALIAS espresso_pyx_flags)
target_compile_options(
  espresso_pyx_flags
  INTERFACE
    $<$<NOT:$<CXX_COMPILER_ID:Intel,IntelLLVM>>:-Wno-pedantic>
    $<$<NOT:$<CXX_COMPILER_ID:Intel,IntelLLVM>>:-Wno-deprecated-declarations>
    $<$<CXX_COMPILER_ID:Intel,IntelLLVM>:-diag-disable=1224>
    $<$<CXX_COMPILER_ID:GNU>:-Wno-cpp>
    $<$<CXX_COMPILER_ID:GNU>:-Wno-strict-aliasing>
    $<$<CXX_COMPILER_ID:GNU>:-Wno-maybe-uninitialized>
    $<$<CXX_COMPILER_ID:GNU>:-Wno-volatile>
    $<$<CXX_COMPILER_ID:Clang,AppleClang,IntelLLVM>:-Wno-sometimes-uninitialized>
    $<$<CXX_COMPILER_ID:Clang,AppleClang,IntelLLVM>:-Wno-\#warnings>
    -Wno-unused-variable)

# Configure, compile and install Cython files
foreach(cython_file ${cython_SRC})
  get_filename_component(basename ${cython_file} NAME_WE)
  file(RELATIVE_PATH relpath ${CMAKE_CURRENT_SOURCE_DIR} ${cython_file})
  if(basename STREQUAL "code_info")
    file(RELATIVE_PATH relpath ${CMAKE_CURRENT_BINARY_DIR} ${cython_file})
  endif()
  get_filename_component(relpath ${relpath} DIRECTORY)
  if(relpath STREQUAL "")
    string(CONCAT outputpath ${CMAKE_CURRENT_BINARY_DIR} "/" ${basename} ".cpp")
  else()
    string(CONCAT outputpath ${CMAKE_CURRENT_BINARY_DIR} "/" ${relpath} "/"
                  ${basename} ".cpp")
  endif()
  if(basename STREQUAL "")
    message(FATAL_ERROR "Internal error empty basename of file ${cython_file}")
  else()
    add_custom_command(
      OUTPUT ${outputpath}
      COMMAND
        ${CYTHON_EXECUTABLE}
        $<$<BOOL:${ESPRESSO_WARNINGS_ARE_ERRORS}>:--warning-errors> -3 --cplus
        --directive embedsignature=True --directive binding=True -I
        ${CMAKE_CURRENT_SOURCE_DIR} -I ${CMAKE_CURRENT_BINARY_DIR}
        ${cython_file} -o ${outputpath}
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..
      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/myconfig.pxi ${cython_file}
              ${cython_HEADER})
    set(target "espressomd_${basename}")
    add_library(${target} SHARED ${outputpath})
    if(NOT "${relpath}" STREQUAL "")
      set_target_properties(${target} PROPERTIES PREFIX "${relpath}/")
    else()
      set_target_properties(${target} PROPERTIES PREFIX "")
    endif()
    set_target_properties(${target} PROPERTIES OUTPUT_NAME ${basename})
    if(APPLE)
      set_target_properties(
        ${target} PROPERTIES SUFFIX ".so" LINK_FLAGS
                                          "-undefined dynamic_lookup")
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
    install(TARGETS ${target}
            LIBRARY DESTINATION ${ESPRESSO_INSTALL_PYTHON}/espressomd)
  endif()
endforeach()

target_link_libraries(espressomd_profiler PRIVATE espresso::profiler)

# Configure Python files
foreach(auxfile ${cython_AUX})
  get_filename_component(filename ${auxfile} NAME)
  file(RELATIVE_PATH relpath ${CMAKE_CURRENT_SOURCE_DIR} ${auxfile})
  get_filename_component(relpath ${relpath} DIRECTORY)
  string(CONCAT outputpath ${CMAKE_CURRENT_BINARY_DIR} "/" ${relpath} "/"
                ${filename})
  add_custom_command(TARGET espressomd COMMAND ${CMAKE_COMMAND} -E copy
                                               ${auxfile} ${outputpath})
endforeach(auxfile)

# Install Python files
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        DESTINATION ${ESPRESSO_INSTALL_PYTHON} FILES_MATCHING PATTERN "*.py"
        PATTERN "CMakeFiles" EXCLUDE)
