#
# Copyright (C) 2012,2013,2014,2015,2016 The ESPResSo project
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
proc writevtk {filename {types "all"}} {
	set max_pid [setmd max_part]
	set n 0
	set fp [open $filename "w"]

	for { set pid 0 } { $pid <= $max_pid } { incr pid } {
		foreach type $types {
			if {[part $pid print type] == $type || ([part $pid print type] != "na" && $type == "all")} then {
				incr n
			}
		}
	}

	puts $fp "# vtk DataFile Version 2.0\nparticles\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS $n floats"

	for { set pid 0 } { $pid <= $max_pid } { incr pid } {
		foreach type $types {
			if {[part $pid print type] == $type || ([part $pid print type] != "na" && $type == "all")} then {
				set xpos [expr [lindex [part $pid print folded_pos] 0]]
				set ypos [expr [lindex [part $pid print folded_pos] 1]]
				set zpos [expr [lindex [part $pid print folded_pos] 2]]
				puts $fp "$xpos $ypos $zpos"
			}
		}
	}

	puts $fp "POINT_DATA $n"
	puts $fp "SCALARS velocity float 3"
	puts $fp "LOOKUP_TABLE default"

	for { set pid 0 } { $pid <= $max_pid } { incr pid } {
		foreach type $types {
			if {[part $pid print type] == $type || ([part $pid print type] != "na" && $type == "all")} then {
				set xvel [expr [lindex [part $pid print v] 0]]
				set yvel [expr [lindex [part $pid print v] 1]]
				set zvel [expr [lindex [part $pid print v] 2]]
				puts $fp "$xvel $yvel $zvel"
			}
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

#---------------------------------------------------------------------------------------------------------------
# writes .vtk for rhomboid obstacle
proc output_vtk_rhomboid { args } {

    # expected parameters: 
    # corX, corY, corZ - corner
    # aX, aY, aZ - length vector 
    # bX, bY, bZ - width vector 
    # cX, cY, cZ - height vector
    # rhomFile - file to write into

    set corX [lindex $args 0]
    set corY [lindex $args 1]
    set corZ [lindex $args 2]
    set aX [lindex $args 3]
    set aY [lindex $args 4]
    set aZ [lindex $args 5]
    set bX [lindex $args 6]
    set bY [lindex $args 7]
    set bZ [lindex $args 8]
    set cX [lindex $args 9]
    set cY [lindex $args 10]
    set cZ [lindex $args 11]
    set rhomFile [lindex $args 12]

    # write output
    set f [open $rhomFile "w"]

    puts $f "# vtk DataFile Version 3.0"
    puts $f "Data"
    puts $f "ASCII"
    puts $f "DATASET POLYDATA"
    puts $f "POINTS 8 float"

    puts $f "$corX $corY $corZ"
    puts $f "[expr $corX + $aX] [expr $corY + $aY] [expr $corZ + $aZ]"
    puts $f "[expr $corX + $aX + $bX] [expr $corY + $aY + $bY] [expr $corZ + $aZ + $bZ]"
    puts $f "[expr $corX + $bX] [expr $corY + $bY] [expr $corZ + $bZ]"

    puts $f "[expr $corX + $cX] [expr $corY + $cY] [expr $corZ + $cZ]"
    puts $f "[expr $corX + $aX + $cX] [expr $corY + $aY + $cY] [expr $corZ + $aZ + $cZ]"
    puts $f "[expr $corX + $aX + $bX + $cX] [expr $corY + $aY + $bY + $cY] [expr $corZ + $aZ + $bZ + $cZ]"
    puts $f "[expr $corX + $bX + $cX] [expr $corY + $bY + $cY] [expr $corZ + $bZ + $cZ]"

    puts $f "POLYGONS 6 30"
    puts $f "4 0 1 2 3"
    puts $f "4 4 5 6 7"
    puts $f "4 0 1 5 4"
    puts $f "4 2 3 7 6"
    puts $f "4 0 4 7 3"
    puts $f "4 1 2 6 5"

    close $f
}
#---------------------------------------------------------------------------------------------------------------
# vrites .vtk for cylinder obstacle
proc output_vtk_cylinder { args } {

    # expected parameters: 
    # cX, cY, cZ - center
    # nX, nY, nZ - normal vector, for now only [0,0,1] 
    # L - half cylinder length 
    # r - radius
    # n - number of faces on the circumference (the higher the number, the smoother the cylinder)
    # cylFile - file to write into

    set cX [lindex $args 0]
    set cY [lindex $args 1]
    set cZ [lindex $args 2]
    set nX [lindex $args 3]
    set nY [lindex $args 4]
    set nZ [lindex $args 5]
    set r [lindex $args 6]
    set L [lindex $args 7]
    set n [lindex $args 8]
    set cylFile [lindex $args 9]

    # write output
    set f [open $cylFile "w"]

    set check_normal 1
    if { $nX != 0.0 } { set check_normal 0 }
    if { $nY != 0.0 } { set check_normal 0 }
    if { $nZ == 0.0 } { set check_normal 0 }
    if { $check_normal == 0 } { 
	puts "This type of cylinder is not supported yet." 
    } else {
	if { $nZ != 1.0 } { set nZ 1.0 }
  
	# set points on the circumference
	set pi 3.14159265359
	set alpha [expr 2*$pi/$n]
	set points [expr 2*$n]

	# get center P1 of bottom circle
	set p1X [expr $cX-$L*$nX]
	set p1Y [expr $cY-$L*$nY]
	set p1Z [expr $cZ-$L*$nZ]

	puts $f "# vtk DataFile Version 3.0"
	puts $f "Data"
	puts $f "ASCII"
	puts $f "DATASET POLYDATA"
	puts $f "POINTS $points float"

	for {set i 0} {$i < $n} {incr i} {
	    puts $f "[expr $p1X+$r*cos($i*$alpha)] [expr $p1Y+$r*sin($i*$alpha)] $p1Z"
	}

	for {set i 0} {$i < $n} {incr i} {
	    puts $f "[expr $p1X+$r*cos($i*$alpha)] [expr $p1Y+$r*sin($i*$alpha)] [expr $p1Z+2*$L*$nZ]"
	}

	puts $f "POLYGONS [expr $n+2] [expr 5*$n+($n+1)*2]"

	# writing the bottom "circle"
	puts -nonewline $f "$n "
	for {set i 0} {$i < [expr $n-1]} {incr i} {
	    puts -nonewline $f "$i "
	}
	puts $f "[expr $n-1]"

	# writing the top "circle"
	puts -nonewline $f "$n "
	for {set i 0} {$i < [expr $n-1]} {incr i} {
	    puts -nonewline $f "[expr $i+$n] "
	}
	puts $f "[expr 2*$n-1]"

	# writing the side rectangles
	for {set i 0} {$i < [expr $n-1]} {incr i} {
	    puts $f "4 $i [expr $i+1] [expr $i+$n+1] [expr $i+$n]"
	}    
	puts $f "4 [expr $n-1] 0 $n [expr 2*$n-1]"

	close $f
    }
}
