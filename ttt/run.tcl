set box_l 40
set sphere_r 10
set temp 0.1
set force 0.005
set gpu 1
set step 100


setmd skin 0.1
setmd box_l $box_l $box_l $box_l
setmd time_step 1.0
cellsystem domain_decomposition -no_verlet_list

if {$gpu == 0} {
	lbfluid den 1 agrid 1 tau 1 visc 1 ext_force $force 0 0
} else {
	lbfluid gpu den 1 agrid 1 tau 1 visc 1 ext_force $force 0 0
}

thermostat lb $temp
lbboundary sphere center [expr $box_l/2] [expr $box_l/2] [expr $box_l/2] radius $sphere_r direction outside

if {$gpu == 0} {
	lbfluid print vtk boundary boundary.vtk ;#remove vtk for gnuplot format
	puts "Wrote boundary files"
}

for {set i 0} {1} {incr i $step} {
	if {$gpu == 0} {
		lbfluid print vtk velocity velocity_$i.vtk ;#remove vtk for gnuplot format
	} else {
		lbprint velocity vtk velocity_$i.vtk
	}
	
  puts -nonewline "integrating $i"
	puts "th step"
	flush stdout
	
  integrate $step
}
