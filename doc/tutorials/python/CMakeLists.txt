# Copyright (C) 2009,2010,2011,2012,2015,2016 The ESPResSo project
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


### Here: add new directory
add_subdirectory(00-building_espresso)
add_subdirectory(01-lennard_jones)
add_subdirectory(02-charged_system)
add_subdirectory(05-raspberry_electrophoresis)
add_subdirectory(08-visualization)


### here: add the appropriate tutorial target after DEPENDS
add_custom_target(tutorials_py DEPENDS tutorials_py_00 tutorials_py_01 tutorials_py_02 tutorials_py_05 tutorials_py_08)
#                                       tutorials_py_03                )
