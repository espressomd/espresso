#
# Copyright (C) 2015-2024 The ESPResSo project
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

find_package(Doxygen)
if(DOXYGEN_FOUND)
  add_custom_command(
    OUTPUT doxy-features
    COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/gen_doxyconfig.py
            ${CMAKE_SOURCE_DIR}/src/config doxy-features
    DEPENDS espresso::config)

  set(DOXYGEN_BIB_IN ${CMAKE_SOURCE_DIR}/doc/bibliography.bib)
  set(DOXYGEN_BIB_OUT ${CMAKE_CURRENT_BINARY_DIR}/bibliography.bib)

  add_custom_command(
    OUTPUT ${DOXYGEN_BIB_OUT}
    COMMAND
      ${CMAKE_COMMAND} -DINPUT=${DOXYGEN_BIB_IN} -DOUTPUT=${DOXYGEN_BIB_OUT} -P
      ${CMAKE_CURRENT_SOURCE_DIR}/bibtex_preprocess.cmake
    DEPENDS ${DOXYGEN_BIB_IN})

  set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
  set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
  set(DOXYGEN_CITELIST ${CMAKE_CURRENT_BINARY_DIR}/html/citelist.html)

  configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

  add_custom_target(
    doxygen COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
    COMMAND ${CMAKE_COMMAND} -DBIBLIOGRAPHY=${DOXYGEN_CITELIST} -P
            ${CMAKE_CURRENT_SOURCE_DIR}/bibtex_postprocess.cmake
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen"
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doxy-features ${DOXYGEN_BIB_OUT}
    VERBATIM)
endif(DOXYGEN_FOUND)
