#
# Copyright (C) 2009-2022 The ESPResSo project
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
#

add_custom_target(tutorials)
add_custom_target(tutorials_html)
add_custom_target(tutorials_python)

configure_file(Readme.md ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(convert.py ${CMAKE_CURRENT_BINARY_DIR})

# Convert a tutorial to the Python and HTML formats. Make these files
# dependencies of targets `tutorials_html` and `tutorials_python`.
function(ESPRESSO_ADD_TUTORIAL)
  cmake_parse_arguments(NB_EXPORT "HTML_RUN" "FILE;TARGET;SUFFIX"
                        "VAR_SUBST;ADD_SCRIPTS;DEPENDS" ${ARGN})

  # create target if it doesn't already exist
  set(TUTORIAL_TARGET ${NB_EXPORT_TARGET})
  if(NOT TARGET ${TUTORIAL_TARGET})
    add_custom_target(${TUTORIAL_TARGET})
  endif()

  # deploy dependencies
  set(TARGET_DEPENDENCIES "")
  foreach(RELPATH ${NB_EXPORT_DEPENDS})
    configure_file(${RELPATH} ${RELPATH} COPYONLY)
    list(APPEND TARGET_DEPENDENCIES ${CMAKE_CURRENT_BINARY_DIR}/${RELPATH})
  endforeach()

  set(NB_FILE ${CMAKE_CURRENT_SOURCE_DIR}/${NB_EXPORT_FILE})
  set(IPYNB_FILE ${CMAKE_CURRENT_BINARY_DIR}/${NB_EXPORT_FILE})
  set(DEPENDENCY_OF_TARGET "${NB_EXPORT_TARGET}")
  if(DEFINED NB_EXPORT_SUFFIX AND NOT "${NB_EXPORT_SUFFIX}" STREQUAL "")
    set(NB_EXPORT_TARGET "${NB_EXPORT_TARGET}_${NB_EXPORT_SUFFIX}")
  endif()
  cmake_path(GET NB_FILE STEM NB_FILE_STEM)
  cmake_path(GET NB_FILE EXTENSION NB_FILE_EXT)
  set(HTML_FILE "${CMAKE_CURRENT_BINARY_DIR}/${NB_FILE_STEM}.html")
  set(PY_FILE "${CMAKE_CURRENT_BINARY_DIR}/${NB_FILE_STEM}.py")

  if(NB_EXPORT_HTML_RUN)
    set(NB_FILE_RUN "${NB_FILE_STEM}.run${NB_FILE_EXT}")
    add_custom_command(
      OUTPUT ${NB_FILE_RUN}
      DEPENDS ${NB_FILE} ${NB_EXPORT_ADD_SCRIPTS} ${CMAKE_BINARY_DIR}/pypresso
              ${CMAKE_BINARY_DIR}/doc/tutorials/convert.py
              ${CMAKE_SOURCE_DIR}/maintainer/parsing/importlib_wrapper.py
      COMMAND
        ${CMAKE_BINARY_DIR}/pypresso
        ${CMAKE_BINARY_DIR}/doc/tutorials/convert.py ci --execute
        --prepare-for-html --input ${NB_FILE} --output ${NB_FILE_RUN}
        --substitutions ${NB_EXPORT_VAR_SUBST} --scripts
        ${NB_EXPORT_ADD_SCRIPTS})
  else()
    set(NB_FILE_RUN ${NB_FILE})
  endif()

  add_custom_command(
    OUTPUT ${IPYNB_FILE}
    DEPENDS ${NB_FILE} ${CMAKE_BINARY_DIR}/pypresso
            ${CMAKE_BINARY_DIR}/doc/tutorials/convert.py
    COMMAND ${CMAKE_COMMAND} -E copy ${NB_FILE} ${IPYNB_FILE}
    COMMAND
      ${CMAKE_BINARY_DIR}/pypresso ${CMAKE_BINARY_DIR}/doc/tutorials/convert.py
      cells --to-md ${IPYNB_FILE})

  add_custom_command(
    OUTPUT ${HTML_FILE} DEPENDS ${NB_FILE_RUN} ${NB_EXPORT_ADD_SCRIPTS}
    COMMAND ${IPYTHON_EXECUTABLE} nbconvert --to "html" --output ${HTML_FILE}
            ${NB_FILE_RUN})

  add_custom_command(
    OUTPUT ${PY_FILE} DEPENDS ${NB_FILE}
    COMMAND ${IPYTHON_EXECUTABLE} nbconvert --to "python" --output ${PY_FILE}
            ${NB_FILE})

  add_custom_target("${NB_EXPORT_TARGET}_deps" DEPENDS ${TARGET_DEPENDENCIES}
                                                       ${IPYNB_FILE})
  add_custom_target("${NB_EXPORT_TARGET}_html" DEPENDS ${HTML_FILE}
                                                       ${TUTORIAL_TARGET})
  add_custom_target("${NB_EXPORT_TARGET}_python" DEPENDS ${PY_FILE}
                                                         ${TUTORIAL_TARGET})
  add_dependencies(${TUTORIAL_TARGET} "${NB_EXPORT_TARGET}_deps")
  add_dependencies(tutorials ${TUTORIAL_TARGET})
  add_dependencies(tutorials_html "${NB_EXPORT_TARGET}_html")
  add_dependencies(tutorials_python "${NB_EXPORT_TARGET}_python")
endfunction()

# Here: add new directory
add_subdirectory(lennard_jones)
add_subdirectory(error_analysis)
add_subdirectory(langevin_dynamics)
add_subdirectory(charged_system)
add_subdirectory(polymers)
add_subdirectory(lattice_boltzmann)
add_subdirectory(raspberry_electrophoresis)
add_subdirectory(active_matter)
add_subdirectory(electrokinetics)
add_subdirectory(visualization)
add_subdirectory(ferrofluid)
add_subdirectory(constant_pH)
add_subdirectory(widom_insertion)
add_subdirectory(electrodes)
add_subdirectory(grand_canonical_monte_carlo)
add_subdirectory(mlip)
add_subdirectory(boltzmann_inversion)
