set box_l 40
set temp 0
set force 0.002
set nue 1
set gpu 0
set step 10


setmd skin 0.1
setmd box_l $box_l $box_l $box_l
setmd time_step 1.0
cellsystem domain_decomposition -no_verlet_list

if {$gpu == 0} {
	lbfluid den 1 agrid 1 tau 1 visc $nue ext_force 0 0 $force
} else {
	lbfluid gpu den 1 agrid 1 tau 1 visc $nue ext_force 0 0 $force
}

thermostat lb $temp

constraint rhomboid corner 20 20 20 a 5 0 0 b 0 5 0 c 0 0 5 direction outside type 1

if {$gpu == 0} {
	lbboundary sphere center 20 20 20 radius 10 direction outside
	lbfluid print vtk boundary boundary.vtk ;#remove vtk for gnuplot format
	puts "Wrote boundary file"
}

for {set i 0} {1} {incr i $step} {
	if {$gpu == 0} {
		lbfluid print vtk velocity velocity_$i.vtk ;#remove vtk for gnuplot format
	} else {
		lbprint velocity vtk velocity_$i.vtk
	}
	
	if {$i % 100 == 0} {
		puts -nonewline "integrating $i"
		puts "th step"
		flush stdout
	}
	
  integrate $step
}
