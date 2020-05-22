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

# Like `option()`, but takes an additional list of valid values that is used
# to check the user input, generate a help string with all options, and set
# the list of possible values in ccmake/cmake-gui (e.g. in ccmake, one can
# cycle through the list of values with [enter] instead of typing them).
macro(option_enum)
  cmake_parse_arguments(option_enum "" "varname;help_text;default_value" "possible_values" ${ARGN})
  # Process user input
  set(option_enum_input_value "${${option_enum_varname}}")
  if(NOT option_enum_input_value IN_LIST option_enum_possible_values)
    if(option_enum_input_value)
      message(WARNING "Unknown ${option_enum_help_text} '${option_enum_input_value}', defaulting to '${option_enum_default_value}'")
    endif()
    set(option_enum_input_value ${option_enum_default_value})
  endif()
  # Declare variable with a help string
  string(REPLACE ";" " " option_enum_possible_values_concat "${option_enum_possible_values}")
  set(${option_enum_varname} ${option_enum_input_value} CACHE STRING "Choose the ${option_enum_help_text}, options are: ${option_enum_possible_values_concat}" FORCE)
  # Set the possible values for ccmake and cmake-gui
  set_property(CACHE ${option_enum_varname} PROPERTY STRINGS ${option_enum_possible_values})
endmacro()
