#
# Copyright (C) 2016-2022 The ESPResSo project
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

include(FindPackageHandleStandardArgs)

set(SPHINX_EXECUTABLE ${Python_EXECUTABLE} -m sphinx)

execute_process(
  COMMAND ${SPHINX_EXECUTABLE} --version OUTPUT_VARIABLE QUERY_VERSION_OUT
  ERROR_VARIABLE QUERY_VERSION_ERR RESULT_VARIABLE QUERY_VERSION_RESULT)

if(NOT QUERY_VERSION_RESULT)
  # Sphinx switched at some point from returning their version on stdout to
  # printing it at stderr. Since we do not know their version yet, we use stdout
  # if it matches a version regex, or stderr otherwise.
  if(QUERY_VERSION_OUT MATCHES "[0-9]+\.[0-9.]+")
    set(QUERY_VERSION "${QUERY_VERSION_OUT}")
  else()
    set(QUERY_VERSION "${QUERY_VERSION_ERR}")
  endif()

  string(REGEX MATCH "[0-9]+\.[0-9.]+" SPHINX_VERSION "${QUERY_VERSION}")

  if("${SPHINX_VERSION}" VERSION_LESS "1.7")
    set(SPHINX_API_DOC_EXE ${Python_EXECUTABLE} -m sphinx.apidoc)
  else()
    set(SPHINX_API_DOC_EXE ${Python_EXECUTABLE} -m sphinx.ext.apidoc)
  endif()
endif()

find_package_handle_standard_args(
  Sphinx REQUIRED_VARS SPHINX_EXECUTABLE SPHINX_API_DOC_EXE
  VERSION_VAR SPHINX_VERSION)

mark_as_advanced(SPHINX_EXECUTABLE)
mark_as_advanced(SPHINX_API_DOC_EXE)
