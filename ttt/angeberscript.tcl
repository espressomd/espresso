set box_x 160
set box_y 40
set box_z 10
set temp 3.0
set N 2000


setmd skin 0.1
setmd box_l $box_x $box_y $box_z
setmd time_step 0.008
cellsystem domain_decomposition -no_verlet_list

lbfluid gpu den 1 agrid 1. tau 0.008 visc 1.0 friction 1. ext_force 1 0 0
thermostat lb $temp

constraint rhomboid corner 5 10 0 a [expr $box_y/4] 0 0 b 0 [expr $box_y/2] 0 c 0 0 $box_z direction outside type 0
lbboundary rhomboid corner 5 10 0 a [expr $box_y/4] 0 0 b 0 [expr $box_y/2] 0 c 0 0 $box_z direction outside

constraint rhomboid corner -1 -1 -1 a [expr $box_x+2] 0 0 b 0 2 0 c 0 0 [expr $box_z+2] direction outside type 0
lbboundary rhomboid corner -1 -1 -1 a [expr $box_x+2] 0 0 b 0 2 0 c 0 0 [expr $box_z+2] direction outside

constraint rhomboid corner -1 [expr $box_y+1] -1 a [expr $box_x+2] 0 0 b 0 -2 0 c 0 0 [expr $box_z+2] direction outside type 0
lbboundary rhomboid corner -1 [expr $box_y+1] -1 a [expr $box_x+2] 0 0 b 0 -2 0 c 0 0 [expr $box_z+2] direction outside

#lbfluid print vtk boundary boundary.vtk

if {$N > 0} {
  set x [expr rand()*$box_x]
  set y [expr rand()*$box_y]
  set z [expr rand()*$box_z]
	    	
  while {[constraint mindist_position $x $y $z] < 1.1} {
	  set x [expr rand()*$box_x]
	  set y [expr rand()*$box_y]
	  set z [expr rand()*$box_z]
  }
		
  part 0 pos $x $y $z type 0
  puts "placed particle 0"
}

for {set id 1} {$id < $N} {incr id} {
	set dist 0.
	
	while { $dist < 1.1} {
	  set x [expr rand()*$box_x]
		set y [expr rand()*$box_y]
		set z [expr rand()*$box_z]
	  	
	  while {[constraint mindist_position $x $y $z] < 1.1} {
		  set x [expr rand()*$box_x]
		  set y [expr rand()*$box_y]
		  set z [expr rand()*$box_z]
		}
		
		part $id pos $x $y $z type 0
		set dist [analyze mindist]
	}
	
  puts "placed particle $id"
}

inter 0 0 lennard-jones 1. 1. 1.1225 0.25 0

prepare_vmd_connection "vmd" 2000

for {set i 0} {1} {incr i} {
	if {$i % 10 == 0} {
		#lbfluid print vtk velocity velocity_$i.vtk	
		puts -nonewline "integrating $i"
		puts "th step"
		flush stdout
	}
	
  integrate 10
  imd positions
}

#lbfluid print vtk velocity velocity.vtk
