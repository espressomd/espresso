#
# Copyright (C) 2013-2018 The ESPResSo project
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

add_custom_command(OUTPUT gen_pxiconfig.cpp
                   COMMAND ${PYTHON_EXECUTABLE}
                   ${CMAKE_CURRENT_SOURCE_DIR}/gen_pxiconfig.py
                   ${CMAKE_SOURCE_DIR}/src/config/features.def ${CMAKE_CURRENT_BINARY_DIR}/gen_pxiconfig.cpp
                   DEPENDS ${CMAKE_SOURCE_DIR}/src/config/features.def
                  )

add_executable(gen_pxiconfig gen_pxiconfig.cpp)
target_link_libraries(gen_pxiconfig EspressoConfig)

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/myconfig.pxi
                   COMMAND ${CMAKE_CURRENT_BINARY_DIR}/gen_pxiconfig > ${CMAKE_CURRENT_BINARY_DIR}/myconfig.pxi
                   DEPENDS gen_pxiconfig
                   )

add_custom_target(espressomd)
add_custom_command(OUTPUT code_info.pyx
        COMMAND ${PYTHON_EXECUTABLE} gen_code_info.py ${CMAKE_SOURCE_DIR}/src/config/features.def ${CMAKE_CURRENT_BINARY_DIR}/code_info.pyx
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/myconfig.pxi
        )

file(GLOB cython_SRC *.pyx)
# Make the cython_SRC, cython_HEADER and cython_AUX a cached variable to
# be able to extend it in the subdirectories.
SET(cython_SRC  "${cython_SRC}" CACHE INTERNAL "cython_SRC")
file(GLOB cython_HEADER *.pxd)
SET(cython_HEADER  "${cython_HEADER}" CACHE INTERNAL "cython_HEADER")
file(GLOB cython_AUX *.py)
configure_file(MDA_ESP/__init__.py ${CMAKE_CURRENT_BINARY_DIR}/MDA_ESP/__init__.py COPYONLY)
SET(cython_AUX  "${cython_AUX}" CACHE INTERNAL "cython_AUX")

add_subdirectory(io)
add_subdirectory(core)

include(CythonModule)

function(add_espresso_cython_module)
  cmake_parse_arguments(CY "" "PYX" "DEPENDS" ${ARGN})
  add_cython_module(TARGET_DIR ${CMAKE_CURRENT_BINARY_DIR} PYX ${CY_PYX} TARGET_PREFIX "espressomd" DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/myconfig.pxi ${cython_HEADER} ${CY_DEPENDENDS} LINK_WITH EspressoConfig EspressoCore EspressoScriptInterface)
endfunction()

add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/actors.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/analyze.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/cellsystem.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/collision_detection.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/cuda_init.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/diamond.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/electrokinetics.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/electrostatic_extensions.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/electrostatics.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/galilei.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/globals.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/_init.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/integrate.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/interactions.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/lb.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/magnetostatic_extensions.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/magnetostatics.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/minimize_energy.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/particle_data.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/polymer.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/profiler.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/reaction_ensemble.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/scafacos.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/script_interface.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/system.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/thermostat.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/utils.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/version.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/visualization_mayavi.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/visualization_opengl.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_SOURCE_DIR}/visualization.pyx)
add_espresso_cython_module(PYX ${CMAKE_CURRENT_BINARY_DIR}/code_info.pyx)


target_link_libraries(espressomd_profiler PRIVATE Profiler)

foreach(auxfile ${cython_AUX})
  get_filename_component(filename ${auxfile} NAME)
  file(RELATIVE_PATH relpath ${CMAKE_CURRENT_SOURCE_DIR} ${auxfile})
  get_filename_component(relpath ${relpath} DIRECTORY)
  string(CONCAT outputpath ${CMAKE_CURRENT_BINARY_DIR} "/" ${relpath} "/" ${filename})
  add_custom_command(TARGET espressomd
                     COMMAND ${CMAKE_COMMAND} -E
                       copy ${auxfile} ${outputpath}
                    )
endforeach(auxfile)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} DESTINATION ${PYTHON_INSTDIR}
	PATTERN "CMakeFiles" EXCLUDE PATTERN "*.cpp" EXCLUDE
	PATTERN "*.cmake" EXCLUDE PATTERN "Makefile" EXCLUDE) 

