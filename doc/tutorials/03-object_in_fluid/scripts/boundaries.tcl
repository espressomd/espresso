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

#-----------------------------------------------------------------------------------------------------
# define your own walls and boundaries here
#
# remember that 
# output_vtk_* writes boundary for visualisation later
# constraint sets up boundary for objects 
# and lbboundary sets up boundary for fluid

# wall - bottom
set corX 0; set corY 0; set corZ 0;
set aX $boxX; set aY 0; set aZ 0;
set bX 0; set bY $boxY; set bZ 0;
set cX 0; set cY 0; set cZ 1;
set rhomFile "output/wallbottom.vtk"
output_vtk_rhomboid $corX $corY $corZ $aX $aY $aZ $bX $bY $bZ $cX $cY $cZ $rhomFile 
constraint rhomboid corner $corX $corY $corZ a $aX $aY $aZ b $bX $bY $bZ c $cX $cY $cZ direction 1 type 10 
lbboundary rhomboid corner $corX $corY $corZ a $aX $aY $aZ b $bX $bY $bZ c $cX $cY $cZ direction 1

# wall - top
set corX 0; set corY 0; set corZ [expr $boxZ - 1];
set aX $boxX; set aY 0; set aZ 0;
set bX 0; set bY $boxY; set bZ 0;
set cX 0; set cY 0; set cZ 1;
set rhomFile "output/walltop.vtk"
output_vtk_rhomboid $corX $corY $corZ $aX $aY $aZ $bX $bY $bZ $cX $cY $cZ $rhomFile 
constraint rhomboid corner $corX $corY $corZ a $aX $aY $aZ b $bX $bY $bZ c $cX $cY $cZ direction 1 type 10 
lbboundary rhomboid corner $corX $corY $corZ a $aX $aY $aZ b $bX $bY $bZ c $cX $cY $cZ direction 1

# wall - front side
set corX 0; set corY 0; set corZ 1;
set aX $boxX; set aY 0; set aZ 0;
set bX 0; set bY 1; set bZ 0;
set cX 0; set cY 0; set cZ [expr $boxZ - 2];
set rhomFile "output/wallfront.vtk"
output_vtk_rhomboid $corX $corY $corZ $aX $aY $aZ $bX $bY $bZ $cX $cY $cZ $rhomFile 
constraint rhomboid corner $corX $corY $corZ a $aX $aY $aZ b $bX $bY $bZ c $cX $cY $cZ direction 1 type 10 
lbboundary rhomboid corner $corX $corY $corZ a $aX $aY $aZ b $bX $bY $bZ c $cX $cY $cZ direction 1

# wall - back side
set corX 0; set corY [expr $boxY - 1]; set corZ 1;
set aX $boxX; set aY 0; set aZ 0;
set bX 0; set bY 1; set bZ 0;
set cX 0; set cY 0; set cZ [expr $boxZ - 2];
set rhomFile "output/wallback.vtk"
output_vtk_rhomboid $corX $corY $corZ $aX $aY $aZ $bX $bY $bZ $cX $cY $cZ $rhomFile 
constraint rhomboid corner $corX $corY $corZ a $aX $aY $aZ b $bX $bY $bZ c $cX $cY $cZ direction 1 type 10 
lbboundary rhomboid corner $corX $corY $corZ a $aX $aY $aZ b $bX $bY $bZ c $cX $cY $cZ direction 1


# obstacle rhomboid1
set corX 10; set corY 1; set corZ 1;
set aX 8; set aY 0; set aZ 0;
set bX 0; set bY 4; set bZ 0;
set cX 0; set cY 0; set cZ 18;
set rhomFile "output/rhomboid1.vtk"
output_vtk_rhomboid $corX $corY $corZ $aX $aY $aZ $bX $bY $bZ $cX $cY $cZ $rhomFile 
constraint rhomboid corner $corX $corY $corZ a $aX $aY $aZ b $bX $bY $bZ c $cX $cY $cZ direction 1 type 10 
lbboundary rhomboid corner $corX $corY $corZ a $aX $aY $aZ b $bX $bY $bZ c $cX $cY $cZ direction 1

# obstacle cylinder1 
set cX 16; set cY 17; set cZ 10;
set nX 0; set nY 0; set nZ 1;
set L 9
set r 3
set cylFile "output/cylinder1.vtk"
set n 20
output_vtk_cylinder $cX $cY $cZ $nX $nY $nZ $r $L $n $cylFile
constraint cylinder center $cX $cY $cZ axis $nX $nY $nZ radius $r length $L direction 1 type 10 
lbboundary cylinder center $cX $cY $cZ axis $nX $nY $nZ radius $r length $L direction 1

# obstacle rhomboid2
set corX 25; set corY 1; set corZ 1;
set aX 5; set aY 0; set aZ 0;
set bX 0; set bY 20; set bZ 0;
set cX 0; set cY 0; set cZ 10;
set rhomFile "output/rhomboid2.vtk"
output_vtk_rhomboid $corX $corY $corZ $aX $aY $aZ $bX $bY $bZ $cX $cY $cZ $rhomFile 
constraint rhomboid corner $corX $corY $corZ a $aX $aY $aZ b $bX $bY $bZ c $cX $cY $cZ direction 1 type 10 
lbboundary rhomboid corner $corX $corY $corZ a $aX $aY $aZ b $bX $bY $bZ c $cX $cY $cZ direction 1

# obstacle cylinder2 
set cX 37; set cY 10; set cZ 10;
set nX 0; set nY 0; set nZ 1;
set L 9
set r 3
set cylFile "output/cylinder2.vtk"
set n 20
output_vtk_cylinder $cX $cY $cZ $nX $nY $nZ $r $L $n $cylFile
constraint cylinder center $cX $cY $cZ axis $nX $nY $nZ radius $r length $L direction 1 type 10 
lbboundary cylinder center $cX $cY $cZ axis $nX $nY $nZ radius $r length $L direction 1
