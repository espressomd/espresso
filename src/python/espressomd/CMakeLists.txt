#
# Copyright (C) 2013-2019 The ESPResSo project
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
    ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/gen_pxiconfig.py
    ${CMAKE_SOURCE_DIR}/src/config/features.def
    ${CMAKE_CURRENT_BINARY_DIR}/gen_pxiconfig.cpp
  DEPENDS ${CMAKE_SOURCE_DIR}/src/config/features.def)

add_executable(gen_pxiconfig gen_pxiconfig.cpp)
target_link_libraries(gen_pxiconfig EspressoConfig)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/myconfig.pxi
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/gen_pxiconfig >
          ${CMAKE_CURRENT_BINARY_DIR}/myconfig.pxi DEPENDS gen_pxiconfig)

add_custom_target(espressomd)
add_custom_command(
  OUTPUT code_info.pyx
  COMMAND
    ${PYTHON_EXECUTABLE} gen_code_info.py
    ${CMAKE_SOURCE_DIR}/src/config/features.def
    ${CMAKE_CURRENT_BINARY_DIR}/code_info.pyx
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/myconfig.pxi)

# Make the cython_SRC, cython_HEADER and cython_AUX a cached variable to be able
# to extend it in the subdirectories.
file(GLOB cython_SRC *.pyx)
set(cython_SRC "${cython_SRC}" CACHE INTERNAL "cython_SRC")
file(GLOB cython_HEADER *.pxd)
set(cython_HEADER "${cython_HEADER}" CACHE INTERNAL "cython_HEADER")
file(GLOB cython_AUX *.py)
configure_file(MDA_ESP/__init__.py
               ${CMAKE_CURRENT_BINARY_DIR}/MDA_ESP/__init__.py COPYONLY)
set(cython_AUX "${cython_AUX}" CACHE INTERNAL "cython_AUX")

add_subdirectory(io)

list(APPEND cython_SRC ${CMAKE_CURRENT_BINARY_DIR}/code_info.pyx)
list(REMOVE_DUPLICATES cython_SRC)

add_library(pyx_interface INTERFACE)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  target_compile_options(
    pyx_interface
    INTERFACE -Wno-pedantic -Wno-cpp -Wno-strict-aliasing
              -Wno-maybe-uninitialized -Wno-unused-variable
              -Wno-deprecated-declarations)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL
                                                 "AppleClang")
  target_compile_options(
    pyx_interface
    INTERFACE -Wno-pedantic -Wno-\#warnings -Wno-sometimes-uninitialized
              -Wno-unused-variable -Wno-deprecated-declarations)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  target_compile_options(pyx_interface INTERFACE -wd1224)
else()
  target_compile_options(pyx_interface INTERFACE -Wno-pedantic
                                                 -Wno-unused-variable)
endif()
target_link_libraries(pyx_interface INTERFACE cxx_interface)

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
        ${CYTHON_EXECUTABLE} $<$<BOOL:${WARNINGS_ARE_ERRORS}>:--warning-errors>
        -3 --cplus --directive embedsignature=True --directive binding=True -I
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
    target_link_libraries(${target} PRIVATE EspressoConfig EspressoCore
                                            ScriptInterface)
    target_link_libraries(${target} PRIVATE pyx_interface)
    target_include_directories(${target} SYSTEM PRIVATE ${PYTHON_INCLUDE_DIRS}
                                                        ${NUMPY_INCLUDE_DIR})
    add_dependencies(espressomd ${target})
    install(TARGETS ${target} LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)
  endif()
endforeach()

target_link_libraries(espressomd_profiler PRIVATE Profiler)

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
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} DESTINATION ${PYTHON_INSTDIR}
        FILES_MATCHING PATTERN "*.py" PATTERN "CMakeFiles" EXCLUDE)
