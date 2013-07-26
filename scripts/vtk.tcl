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

#dumps particle positions into a file so that paraview can visualize them
proc writevectorvtk {filename {type "all"}} {
	set max_pid [setmd max_part]
	set nobj 0
	set nbonds 0
	set fp [open $filename "w"]
	set hash {}
## the hash table is structured as follows
## <hash_id> <pid> <x> <y> <z> <has_id of the bonded particle>
##
        set boxhx [expr [lindex [setmd box_l] 0 ] /2.] 
        set boxhy [expr [lindex [setmd box_l] 1 ] /2.] 
        set boxhz [expr [lindex [setmd box_l] 2 ] /2.] 
	for { set pid 0 } { $pid <= $max_pid } { incr pid } {
		if {[part $pid print type] == $type || ([part $pid print type] != "na" && $type == "all")} then {
                   set bondid -1
		   set xpos [expr [lindex [part $pid print folded_pos] 0]]
		   set ypos [expr [lindex [part $pid print folded_pos] 1]]
		   set zpos [expr [lindex [part $pid print folded_pos] 2]]
		   if { [llength [lindex [part $pid print bond ] 0 ] ] > 0  } { 
		       # i.e. there is a bond
                       set bondpid [lindex [lindex [part $pid print bond] 0 0] 1 ]
		       set x2pos [expr [lindex [part $bondpid  print folded_pos] 0]]
		       set y2pos [expr [lindex [part $bondpid  print folded_pos] 1]]
		       set z2pos [expr [lindex [part $bondpid  print folded_pos] 2]]
		       if  { [expr abs($xpos - $x2pos) ] < $boxhx && [expr abs($ypos - $y2pos) ] < $boxhy && [expr abs($zpos - $z2pos) ] < $boxhz } { 
                       # i.e. it doesn't cross pbc, then we represent it
		          	incr nbonds
                          	foreach element $hash  { 
		          		if { [lindex $element 1] == $bondpid } { 
		          			set bondid [lindex $element 0]
		          			break
                                  	}
                       	        }
		       }
		   }
                   lappend hash [list $nobj $pid $bondid ] ; #here $bondid==-1 if there's no bonded particle 
		   incr nobj
		}
	}

	puts $fp "# vtk DataFile Version 2.0\nparticles\nASCII\nDATASET POLYDATA\nPOINTS $nobj floats"

        foreach element $hash { 
			set xpos [expr [lindex [part [lindex $element 1 ] print folded_pos] 0]]
			set ypos [expr [lindex [part [lindex $element 1 ] print folded_pos] 1]]
			set zpos [expr [lindex [part [lindex $element 1 ] print folded_pos] 2]]
			puts $fp "$xpos $ypos $zpos"
        }
        puts $fp "\nVERTICES $nobj [expr int(2*$nobj)]"
	for { set id 0 } { $id < $nobj } { incr id } {  
			puts $fp "1 $id"
        } 
        if { $nbonds > 0 } { 
		puts $fp "\nLINES $nbonds [expr int(3*$nbonds)]"
                foreach element $hash { 
			if { [lindex $element 2] >= 0 } { 
			         puts $fp "2 [lindex $element 0 ] [lindex $element 2 ]"
                        }
                }
	}
        puts $fp "\nPOINT_DATA $nobj"
        puts $fp "VECTORS vector float"
	foreach element $hash { 
		set xv [expr [lindex [part [lindex $element  1 ] print dip] 0]]
		set yv [expr [lindex [part [lindex $element  1 ] print dip] 1]]
		set zv [expr [lindex [part [lindex $element  1 ] print dip] 2]]
		puts $fp "$xv $yv $zv"
	}

	close $fp
}
