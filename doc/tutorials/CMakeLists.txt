# Copyright (C) 2009-2018 The ESPResSo project
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

function(NB_EXPORT)
  cmake_parse_arguments(NB_EXPORT "" "FILE;TUTORIAL;TYPE" "" ${ARGN})
  get_filename_component(BASE_NAME ${NB_EXPORT_FILE} NAME_WE)

  if(${NB_EXPORT_TYPE} STREQUAL "html")
    set(NB_EXT "html")
  elseif(${NB_EXPORT_TYPE} STREQUAL "python")
    set(NB_EXT "py")
  endif()

  set(TARGET "${NB_EXPORT_TUTORIAL}_${NB_EXPORT_TYPE}")
  set(NB_FILE_BASE "${CMAKE_CURRENT_BINARY_DIR}/${BASE_NAME}")
  set(NB_FILE "${NB_FILE_BASE}.${NB_EXT}")

  add_custom_command(
    OUTPUT ${NB_FILE}
    DEPENDS ${NB_EXPORT_FILE}
    COMMAND
    ${IPYTHON_EXECUTABLE} nbconvert
    --to ${NB_EXPORT_TYPE}
    --output "${NB_FILE_BASE}"
    ${NB_EXPORT_FILE}
    )

  add_custom_target(${TARGET} DEPENDS ${NB_FILE})
endfunction(NB_EXPORT)

function(HTML_EXPORT)
  cmake_parse_arguments(HTML_EXPORT "" "FILE;TUTORIAL" "" ${ARGN})
  nb_export(FILE ${HTML_EXPORT_FILE} TUTORIAL ${HTML_EXPORT_TUTORIAL} TYPE "html")
endfunction(HTML_EXPORT)

function(PYTHON_EXPORT)
  cmake_parse_arguments(PYTHON_EXPORT "" "FILE;TUTORIAL" "" ${ARGN})
  nb_export(FILE ${PYTHON_EXPORT_FILE} TUTORIAL ${PYTHON_EXPORT_TUTORIAL} TYPE "python")
endfunction(PYTHON_EXPORT)

### Here: add new directory
add_subdirectory(01-lennard_jones)
add_subdirectory(02-charged_system)
add_subdirectory(04-lattice_boltzmann)
add_subdirectory(05-raspberry_electrophoresis)
add_subdirectory(06-active_matter)
add_subdirectory(07-electrokinetics)
add_subdirectory(08-visualization)
add_subdirectory(09-swimmer_reactions)
add_subdirectory(10-reaction_ensemble)


### here: add the appropriate tutorial target after DEPENDS
add_custom_target(tutorials
  DEPENDS
  tutorials_01
  tutorials_02
  tutorials_04
  tutorials_05
  tutorials_06
  tutorials_07
  tutorials_08
  tutorials_09
  tutorials_10)

add_custom_target(tutorials_html DEPENDS
  tutorials
  tutorials_01_html
  tutorials_02_1_html
  tutorials_02_2_html
  tutorials_04_1_html
  tutorials_04_2_html
  tutorials_04_3_html
  tutorials_04_4_html
  tutorials_05_html
  tutorials_08_html
  )

add_custom_target(tutorials_python_only DEPENDS
  tutorials_01_python
  tutorials_02_1_python
  tutorials_02_2_python
  tutorials_04_1_python
  tutorials_04_2_python
  tutorials_04_3_python
  tutorials_04_4_python
  tutorials_05_python
  tutorials_08_python
  )

add_custom_target(tutorials_python DEPENDS
  tutorials
  tutorials_python_only
  )

configure_file(Readme.rst ${CMAKE_CURRENT_BINARY_DIR}/Readme.rst COPYONLY)
