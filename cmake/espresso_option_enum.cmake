#
# Copyright (C) 2020-2022 The ESPResSo project
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
macro(espresso_option_enum)
  cmake_parse_arguments(espresso_option_enum "" "varname;help_text;default_value" "possible_values" ${ARGN})
  # Process user input
  set(espresso_option_enum_input_value "${${espresso_option_enum_varname}}")
  if(NOT espresso_option_enum_input_value IN_LIST espresso_option_enum_possible_values)
    if(espresso_option_enum_input_value)
      message(WARNING "Unknown ${espresso_option_enum_help_text} '${espresso_option_enum_input_value}', defaulting to '${espresso_option_enum_default_value}'")
    endif()
    set(espresso_option_enum_input_value ${espresso_option_enum_default_value})
  endif()
  # Declare variable with a help string
  string(REPLACE ";" " " espresso_option_enum_possible_values_concat "${espresso_option_enum_possible_values}")
  set(${espresso_option_enum_varname} ${espresso_option_enum_input_value} CACHE STRING "Choose the ${espresso_option_enum_help_text}, options are: ${espresso_option_enum_possible_values_concat}" FORCE)
  # Set the possible values for ccmake and cmake-gui
  set_property(CACHE ${espresso_option_enum_varname} PROPERTY STRINGS ${espresso_option_enum_possible_values})
endmacro()
