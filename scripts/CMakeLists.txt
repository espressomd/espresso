# Copyright (C) 2009,2010,2011,2012 The ESPResSo project
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
add_custom_target(scripts ALL COMMENT "Copying Tcl scripts to builddir")

# copy all tcl files to bin dir
file(GLOB TCL_SCRIPTS *.tcl)
foreach(FILENAME ${TCL_SCRIPTS})
  add_custom_command(TARGET scripts 
    COMMAND ${CMAKE_COMMAND} -E copy_if_different ${FILENAME} ${CMAKE_CURRENT_BINARY_DIR})
  install(FILES ${FILENAME} DESTINATION share/espresso/scripts)
endforeach(FILENAME)

