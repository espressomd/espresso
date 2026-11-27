#
# Copyright (C) 2026 The ESPResSo project
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

# Find the pystencils runtime headers in the Python environment.
#
# The kernels generated with pystencils include headers from
# ``pystencils_runtime/``, which are only distributed with the pystencils Python
# package. This module asks the Python interpreter for their location.
#
# Result variables: ``Pystencils_FOUND``, ``Pystencils_VERSION``,
# ``Pystencils_INCLUDE_DIR``. Imported target: ``Pystencils::runtime``
# (header-only).

if(NOT Python_EXECUTABLE)
  find_package(Python 3.11 QUIET COMPONENTS Interpreter)
endif()

if(Python_EXECUTABLE)
  execute_process(
    COMMAND
      ${Python_EXECUTABLE} -c
      "import pystencils, pystencils.include as inc;print(pystencils.__version__);print(inc.get_pystencils_include_path())"
    OUTPUT_VARIABLE _pystencils_output
    ERROR_VARIABLE _pystencils_error
    RESULT_VARIABLE _pystencils_result
    OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(_pystencils_result EQUAL 0)
    string(REPLACE "\n" ";" _pystencils_output "${_pystencils_output}")
    list(GET _pystencils_output 0 Pystencils_VERSION)
    list(GET _pystencils_output 1 _pystencils_include_dir)
    if(EXISTS "${_pystencils_include_dir}/pystencils_runtime/generic_cpu.hpp")
      set(Pystencils_INCLUDE_DIR "${_pystencils_include_dir}")
    endif()
  elseif(NOT Pystencils_FIND_QUIETLY)
    message(STATUS "Cannot import pystencils with ${Python_EXECUTABLE}:\n"
                   "${_pystencils_error}")
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  Pystencils REQUIRED_VARS Pystencils_INCLUDE_DIR Python_EXECUTABLE
  VERSION_VAR Pystencils_VERSION HANDLE_VERSION_RANGE)

if(Pystencils_FOUND AND NOT TARGET Pystencils::runtime)
  add_library(Pystencils::runtime INTERFACE IMPORTED)
  set_target_properties(
    Pystencils::runtime PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                   "${Pystencils_INCLUDE_DIR}")
endif()

unset(_pystencils_output)
unset(_pystencils_error)
unset(_pystencils_result)
unset(_pystencils_include_dir)
