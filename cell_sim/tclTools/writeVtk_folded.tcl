# description: all tools to write vtk data to visualize in paraview.
# date: 17.02.14
# author: Vera
# functions:
# 	- writevtk: vtk files for particles
# 	- writevtkCell: for Cells
# 	- writevtkPol: for Polymers


#################################################################
# standard Espresso script to write vtks 			#
#################################################################
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


#################################################################
#								#
# Write vtks for a cell with foldet coordinates.		#
# ----------------------------------------------		#
# (triangulars vanish if they are at both sides of the box. 	#
#################################################################

proc writevtkCell {filename {type "all"}} {

	global boxx
	global boxy
	global boxz

	set max_pid [setmd max_part]
	set n 0
	set fp [open $filename "w"]

	for { set pid 0 } { $pid <= $max_pid } { incr pid } {
		if {[part $pid print type] == $type || ([part $pid print type] != "na" && $type == "all")} then {
			incr n
		}
	}

	puts $fp "# vtk DataFile Version 2.0\n3D Unstructured Grid of Triangles\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS $n floats"
	set nodes {}
	for { set pid 0 } { $pid <= $max_pid } { incr pid } {
		if {[part $pid print type] == $type || ([part $pid print type] != "na" && $type == "all")} then {
			set xpos [expr [lindex [part $pid print folded_pos] 0]]
			set ypos [expr [lindex [part $pid print folded_pos] 1]]
			set zpos [expr [lindex [part $pid print folded_pos] 2]]
			puts $fp "$xpos $ypos $zpos"
			lappend nodes "$xpos $ypos $zpos"
		}
	}
	
	
     #  read the given triangles out of the file
     set in [open "/home/mods/Desktop/espresso/espressoCell/buildNeo/sphere/triangles" r]
     fconfigure $in -buffering line
     gets $in data
     set m 0
     set cells {}
     gets $in data
     #set numTriangles 0
     while {$data != ""} {
	  set part1 [lindex $nodes [lindex $data 0]]
	  set part2 [lindex $nodes [lindex $data 1]]
	  set part3 [lindex $nodes [lindex $data 2]]
	  set dist1x [expr ([lindex $part1 0] - [lindex $part2 0])*([lindex $part1 0] - [lindex $part2 0]) ]
	  set dist1y [expr ([lindex $part1 1] - [lindex $part2 1])*([lindex $part1 1] - [lindex $part2 1]) ]
	  set dist1z [expr ([lindex $part1 2] - [lindex $part2 2])*([lindex $part1 2] - [lindex $part2 2]) ]
	  set dist2x [expr ([lindex $part2 0] - [lindex $part3 0])*([lindex $part2 0] - [lindex $part3 0]) ]
	  set dist2y [expr ([lindex $part2 1] - [lindex $part3 1])*([lindex $part2 1] - [lindex $part3 1]) ]
	  set dist2z [expr ([lindex $part2 2] - [lindex $part3 2])*([lindex $part2 2] - [lindex $part3 2]) ]
	  set dist3x [expr ([lindex $part1 0] - [lindex $part3 0])*([lindex $part1 0] - [lindex $part3 0]) ]
	  set dist3y [expr ([lindex $part1 1] - [lindex $part3 1])*([lindex $part1 1] - [lindex $part3 1]) ]
	  set dist3z [expr ([lindex $part1 2] - [lindex $part3 2])*([lindex $part1 2] - [lindex $part3 2]) ]


	  if {	$dist1x < [expr $boxx*$boxx/4] && $dist1y < [expr $boxy*$boxy/4] && $dist1z < [expr $boxz*$boxz/4] && \
		$dist2x < [expr $boxx*$boxx/4] && $dist2y < [expr $boxy*$boxy/4] && $dist2z < [expr $boxz*$boxz/4] && \
		$dist3x < [expr $boxx*$boxx/4] && $dist3y < [expr $boxy*$boxy/4] && $dist3z < [expr $boxz*$boxz/4]} then {
			lappend cells "3 $data"
			incr m }
          gets $in data
     }
     set o [expr 4*$m]
     puts $fp "CELLS $m $o"
     foreach c $cells {puts $fp "$c"}
     close $in
     
     #sets mesh as Triangles
     puts $fp "CELL_TYPES $m"
     for { set pid 0 } { $pid < $m } { incr pid } {
     puts $fp "5"
     }
      
     
	close $fp
}


#################################################################
#								#
# Write vtks for a polymer with foldet coordinates.		#
# ----------------------------------------------		#
# (lines vanish if they are at both sides of the box. 	#
#################################################################

proc writevtkPol {filename {type "all"}} {
	global nummon
	global numpol
	global boxx
	global boxy
	global boxz
	set max_pid [setmd max_part]
	set n 0
	set fp [open $filename "w"]

	for { set pid 0 } { $pid <= $max_pid } { incr pid } {
		if {[part $pid print type] == $type || ([part $pid print type] != "na" && $type == "all")} then {
			incr n
		}
	}

	puts $fp "# vtk DataFile Version 2.0\n3D Unstructured Grid of Lines\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS $n floats"
	set part {}
	for { set pid 0 } { $pid <= $max_pid } { incr pid } {
		if {[part $pid print type] == $type || ([part $pid print type] != "na" && $type == "all")} then {
			set xpos [expr [lindex [part $pid print folded_pos] 0]]
			set ypos [expr [lindex [part $pid print folded_pos] 1]]
			set zpos [expr [lindex [part $pid print folded_pos] 2]]
			puts $fp "$xpos $ypos $zpos"
			lappend part "$xpos $ypos $zpos"
			#puts "added"
		}
	}
	set cells {}
	set num 0
	for {set j 0} {$j < $numpol} {incr j} {	
		for {set i 0} {$i < [expr $nummon-1]} {incr i} {
			set mon1 [lindex $part [expr $j*$nummon +$i]]
			set mon2 [lindex $part [expr $j*$nummon +$i+1]]
			set distx [expr ( [lindex $mon1 0] - [lindex $mon2 0]) * ([lindex $mon1 0] - [lindex $mon2 0] ) ]
			set disty [expr ( [lindex $mon1 1] - [lindex $mon2 1]) * ([lindex $mon1 1] - [lindex $mon2 1] ) ]
			set distz [expr ( [lindex $mon1 2] - [lindex $mon2 2]) * ([lindex $mon1 2] - [lindex $mon2 2] ) ]
			if { $distx < [expr $boxx*$boxx/4] && $disty < [expr $boxy*$boxy/4] && $distz < [expr $boxz*$boxz/4] } then {
				lappend cells "2 [expr $j*$nummon +$i] [expr $j*$nummon +$i+1]"
				incr num}
		}
	}
	puts $fp "CELLS $num [expr $num*3]"
	foreach c $cells { puts $fp "$c"}	
	puts $fp "CELL_TYPES $num"
	for {set i 0} {$i < $num} {incr i} {
		puts $fp "3"
	}
		
	close $fp
}
