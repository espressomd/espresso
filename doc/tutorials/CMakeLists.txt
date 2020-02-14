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
      DEPENDS "${NB_FILE};${NB_EXPORT_ADD_SCRIPTS};${CMAKE_BINARY_DIR}/doc/tutorials/html_runner.py;${CMAKE_BINARY_DIR}/testsuite/scripts/importlib_wrapper.py"
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
    DEPENDS ${HTML_FILE_DEPENDS};${NB_EXPORT_ADD_SCRIPTS}
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
add_subdirectory(11-ferrofluid)
add_subdirectory(12-constant_pH)


### here: add the appropriate tutorial target after DEPENDS
add_custom_target(tutorials DEPENDS
  tutorial_01
  tutorial_02
  tutorial_04
  tutorial_05
  tutorial_06
  tutorial_07
  tutorial_08
  tutorial_11
  tutorial_12
  )

add_custom_target(tutorials_html DEPENDS
  tutorials
  tutorial_01_html
  tutorial_02_1_html
  tutorial_02_2_html
  tutorial_04_1_html
  tutorial_04_2_html
  tutorial_04_3_html
  tutorial_04_4_html
  tutorial_05_html
  tutorial_06_html
  tutorial_07_html
  tutorial_08_html
  tutorial_11_1_html
  tutorial_11_2_html
  tutorial_11_3_html
  tutorial_12_html
  )

add_custom_target(tutorials_python_only DEPENDS
  tutorial_01_python
  tutorial_02_1_python
  tutorial_02_2_python
  tutorial_04_1_python
  tutorial_04_2_python
  tutorial_04_3_python
  tutorial_04_4_python
  tutorial_05_python
  tutorial_07_python
  tutorial_08_python
  tutorial_11_1_python
  tutorial_11_2_python
  tutorial_11_3_python
  tutorial_12_python
  )

add_custom_target(tutorials_python DEPENDS
  tutorials
  tutorials_python_only
  )

configure_file(Readme.rst ${CMAKE_CURRENT_BINARY_DIR}/Readme.rst COPYONLY)
configure_file(html_runner.py ${CMAKE_CURRENT_BINARY_DIR}/html_runner.py)
