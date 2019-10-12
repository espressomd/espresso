# Copyright (C) 2009-2019 The ESPResSo project
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
  cmake_parse_arguments(NB_EXPORT "HTML_RUN" "FILE;TUTORIAL" "VAR_SUBST;ADD_SCRIPTS" ${ARGN})

  set(NB_FILE "${CMAKE_CURRENT_BINARY_DIR}/${NB_EXPORT_FILE}")
  get_filename_component(NB_FILE_BASE ${NB_FILE} NAME_WE)
  get_filename_component(NB_FILE_EXT ${NB_FILE} EXT)
  set(HTML_FILE "${NB_FILE_BASE}.html")
  set(PY_FILE "${NB_FILE_BASE}.py")

  if(${NB_EXPORT_HTML_RUN})
    set(HTML_FILE_DEPENDS "${NB_FILE_BASE}.run${NB_FILE_EXT}")
    add_custom_command(
      OUTPUT "${HTML_FILE_DEPENDS}"
      DEPENDS "${NB_FILE}"
      COMMAND
      "${CMAKE_BINARY_DIR}/pypresso"
      "${CMAKE_BINARY_DIR}/doc/tutorials/html_runner.py"
      "--input" "${NB_FILE}"
      "--output" "${HTML_FILE_DEPENDS}"
      "--substitutions" ${NB_EXPORT_VAR_SUBST}
      "--scripts" ${NB_EXPORT_ADD_SCRIPTS}
      )
  else()
    set(HTML_FILE_DEPENDS "${NB_FILE}")
  endif()

  add_custom_command(
    OUTPUT ${HTML_FILE}
    DEPENDS ${HTML_FILE_DEPENDS}
    COMMAND
    ${IPYTHON_EXECUTABLE} nbconvert
    --to "html"
    --output ${HTML_FILE}
    ${HTML_FILE_DEPENDS}
    )

  add_custom_command(
    OUTPUT ${PY_FILE}
    DEPENDS ${NB_FILE}
    COMMAND
    ${IPYTHON_EXECUTABLE} nbconvert
    --to "python"
    --output ${PY_FILE}
    ${NB_FILE}
    )

  add_custom_target("${NB_EXPORT_TUTORIAL}_html" DEPENDS ${HTML_FILE})
  add_custom_target("${NB_EXPORT_TUTORIAL}_python" DEPENDS ${PY_FILE})
endfunction(NB_EXPORT)


### Here: add new directory
add_subdirectory(01-lennard_jones)
add_subdirectory(02-charged_system)
add_subdirectory(04-lattice_boltzmann)
add_subdirectory(05-raspberry_electrophoresis)
add_subdirectory(06-active_matter)
add_subdirectory(07-electrokinetics)
add_subdirectory(08-visualization)
add_subdirectory(10-reaction_ensemble)
add_subdirectory(11-ferrofluid)
add_subdirectory(12-constant_pH)


### here: add the appropriate tutorial target after DEPENDS
add_custom_target(tutorials DEPENDS
  tutorials_01
  tutorials_02
  tutorials_04
  tutorials_05
  tutorials_06
  tutorials_07
  tutorials_08
  tutorials_10
  tutorials_11
  tutorials_12
  )

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
  tutorials_07_html
  tutorials_08_html
  tutorials_11_1_html
  tutorials_11_2_html
  tutorials_11_3_html
  tutorials_12_html
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
  tutorials_07_python
  tutorials_08_python
  tutorials_11_1_python
  tutorials_11_2_python
  tutorials_11_3_python
  tutorials_12_python
  )

add_custom_target(tutorials_python DEPENDS
  tutorials
  tutorials_python_only
  )

configure_file(Readme.rst ${CMAKE_CURRENT_BINARY_DIR}/Readme.rst COPYONLY)
configure_file(html_runner.py ${CMAKE_CURRENT_BINARY_DIR}/html_runner.py)
