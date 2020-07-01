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

# Like `option()`, but create an extra boolean variable to store whether the
# option was set to its default value or to a user-provided value. With this
# command, the project can be installed with optional dependencies without
# the need to provide a list of CMake flags. Unavailable dependencies will be
# silently ignored. However, if the user specifically requested an optional
# dependency by passing the corresponding CMake flag, the build system has
# the possibility to throw an error if the dependency is unavailable.
#
# Note that when calling CMake again without clearing the build folder,
# variables from the previous CMake call are loaded in memory. For example,
# if the user passed a value to an `option_if_available()` the first time but
# not the second time, the variable will still be flagged as a user-provided
# value in the second CMake call.
macro(option_if_available varname help_text default_value)
  if(NOT DEFINED ${varname}_IS_DEFAULT_VALUE)
    if("${${varname}}" STREQUAL "")
      set(${varname}_IS_DEFAULT_VALUE TRUE CACHE INTERNAL "does ${varname} contain the default value?")
    else()
      set(${varname}_IS_DEFAULT_VALUE FALSE CACHE INTERNAL "does ${varname} contain the default value?")
    endif()
  endif()
  option(${varname} ${help_text} ${default_value})
endmacro()
