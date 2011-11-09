set box_x 200
set box_y 60
set box_z 10
set temp 3.0
set N 4500
set time_step 0.008
set gpu "on"

setmd skin 0.1
setmd box_l $box_x $box_y $box_z
setmd time_step $time_step

if {$gpu == "on"} {
cellsystem domain_decomposition -no_verlet_list
lbfluid gpu dens 1 agrid 1. tau $time_step visc 1.0 friction 1. ext_force 0.2 0 0
} else {
lbfluid cpu dens 1 agrid 1. tau $time_step visc 1.0 friction 1. ext_force 0.2 0 0
}
thermostat lb $temp

#constraint rhomboid corner 5 10 0 a [expr $box_y/2] [ expr $box_y/2 ] 0 b 0 [expr 0.2*$box_y/4] 0 c 0 0 $box_z direction outside type 0
#lbboundary rhomboid corner 5 10 0 a [expr $box_y/2] [ expr $box_y/2 ] 0 b 0 [expr 0.2*$box_y/4] 0 c 0 0 $box_z direction outside

#constraint rhomboid corner 5 10 0 a 5 5 0 b 0 20 0 c 0 0 10 direction outside type 0
#lbboundary rhomboid corner 5 10 0 a 5 5 0 b 0 20 0 c 0 0 10 direction outside

constraint rhomboid corner 5 [expr $box_y/4] 0 a 10 0 0 b 0 [expr $box_y/2] 0 c 0 0 $box_z direction outside type 0
lbboundary rhomboid corner 5 [expr $box_y/4] 0 a 10 0 0 b 0 [expr $box_y/2] 0 c 0 0 $box_z direction outside

constraint rhomboid corner -1 -1 -1 a [expr $box_x+2] 0 0 b 0 2 0 c 0 0 [expr $box_z+2] direction outside type 0
lbboundary rhomboid corner -1 -1 -1 a [expr $box_x+2] 0 0 b 0 2 0 c 0 0 [expr $box_z+2] direction outside

constraint rhomboid corner -1 [expr $box_y+1] -1 a [expr $box_x+2] 0 0 b 0 -2 0 c 0 0 [expr $box_z+2] direction outside type 0
lbboundary rhomboid corner -1 [expr $box_y+1] -1 a [expr $box_x+2] 0 0 b 0 -2 0 c 0 0 [expr $box_z+2] direction outside

#lbfluid print vtk boundary boundary.vtk

#random initial particle distribution
puts "placing $N particles randomly in the box"

for {set id 0} {$id < $N} {incr id} {
  set x [expr rand()*$box_x]
  set y [expr rand()*$box_y]
  set z [expr rand()*$box_z]

	part $id pos $x $y $z type 0
		  	
  while {[constraint mindist_position $x $y $z] < 1.1} {
    set x [expr rand()*$box_x]
	  set y [expr rand()*$box_y]
	  set z [expr rand()*$box_z]

	  part $id pos $x $y $z type 0
  }
}

prepare_vmd_connection "vmd" 2000 1 "1"

#slowly remove force capping
puts "driving up lennard-jones interaction"
inter 0 0 lennard-jones 1. 1. 1.1225 0.25 0

set lj_forcecap 20
inter ljforcecap $lj_forcecap

while {[analyze mindist] < 1.0 && $lj_forcecap < 200} {
  integrate 40
  incr lj_forcecap 10
  inter ljforcecap $lj_forcecap
  imd positions
}
if {$gpu == "on"} {
puts "Tuning cell length and verlet skin"
tune_cells
}
#go
for {set i 0} {1} {incr i} {
	if {$i % 100 == 0} {
		#lbfluid print vtk velocity velocity_$i.vtk	
		puts -nonewline "integrating $i\r"
		#puts "th step"
		flush stdout
	}
	
  integrate 2
  imd positions
}

#lbfluid print vtk velocity velocity.vtk
