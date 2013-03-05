#
# Copyright (C) 2012,2013 The ESPResSo project
# Copyright (C) 2006,2007,2008,2009,2010,2011 Olaf Lenz
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
#############################################################
#                                                           #
# vtk.tcl                                                   #
# =======                                                   #
#                                                           #
# Functions that allow writing VTK files.                   #
#                                                           #
#############################################################

#dumps particle positions into a file so that paraview can visualize them
proc writevtk {filename {type "all"}} {
	set max_pid [setmd max_part]
	set n 0
	set fp [open $filename "w"]

	for { set pid 0 } { $pid <= $max_pid } { incr pid } {
		if {[part $pid print type] == $type || ([part $pid print type] != "na" && $type == "all")} then {
			incr n
		}
	}

	puts $fp "# vtk DataFile Version 2.0\nparticles\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS $n floats"

	for { set pid 0 } { $pid <= $max_pid } { incr pid } {
		if {[part $pid print type] == $type || ([part $pid print type] != "na" && $type == "all")} then {
			set xpos [expr [lindex [part $pid print folded_pos] 0]]
			set ypos [expr [lindex [part $pid print folded_pos] 1]]
			set zpos [expr [lindex [part $pid print folded_pos] 2]]
			puts $fp "$xpos $ypos $zpos"
		}
	}

	close $fp
}
