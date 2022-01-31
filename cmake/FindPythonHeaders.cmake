# Copyright (C) 2020 The ESPResSo project
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

# find the Python C++ headers
execute_process(
  COMMAND ${PYTHON_EXECUTABLE} -c
          "import sysconfig; print(sysconfig.get_path('include'))"
  OUTPUT_VARIABLE PYTHON_INCLUDE_DIRS OUTPUT_STRIP_TRAILING_WHITESPACE)
# find Python installation directory
if(NOT PYTHON_INSTDIR)
  execute_process(
    COMMAND
      ${PYTHON_EXECUTABLE} -c
      "import sysconfig; print(sysconfig.get_path('purelib', vars={'base': ''}).lstrip('/'))"
    OUTPUT_VARIABLE PYTHON_INSTDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
endif(NOT PYTHON_INSTDIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PythonHeaders REQUIRED_VARS
                                  PYTHON_INCLUDE_DIRS PYTHON_INSTDIR)
