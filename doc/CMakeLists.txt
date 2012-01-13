# Copyright (C) 2012 The ESPResSo project
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
find_package(LATEX)
if (LATEX_COMPILER AND BIBTEX_COMPILER AND MAKEINDEX_COMPILER)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/latexit.sh.in 
    ${CMAKE_CURRENT_BINARY_DIR}/latexit.sh)

  macro(LATEXIT BASENAME SOURCEDIRS)
    add_custom_target(${BASENAME} DEPENDS ${BASENAME}.pdf)

    file(GLOB TEX_FILES *.tex *.bib)
    add_custom_command(OUTPUT ${BASENAME}.pdf
      DEPENDS ${TEX_FILES}
      COMMAND /bin/sh ${CMAKE_BINARY_DIR}/doc/latexit.sh ${SOURCEDIRS} ${BASENAME}
      COMMENT "LaTeXing ${BASENAME}.pdf"
      )
  endmacro(LATEXIT)

  add_custom_target(doc)

  add_subdirectory(ug)
  add_subdirectory(dg)
  add_subdirectory(tutorials)
endif(LATEX_COMPILER AND BIBTEX_COMPILER AND MAKEINDEX_COMPILER)

add_subdirectory(doxygen)
add_dependencies(doc doxygen)
