set box_l 40
set temp 1.0
set N 1000


setmd skin 0.1
setmd box_l $box_l $box_l $box_l
setmd time_step 0.01
cellsystem domain_decomposition -no_verlet_list

lbfluid den 1 agrid 2. tau 0.01 visc 1.0 friction 1. ext_force 0 0 .02
thermostat lb $temp
#thermostat langevin 1 1

constraint rhomboid corner 20 20 20 a 5 0 0 b 0 5 0 c 0 0 5 direction outside type 0
lbboundary rhomboid corner 20 20 20 a 5 0 0 b 0 5 0 c 0 0 5 direction outside

lbfluid print vtk boundary boundary.vtk

part 0 pos 40 40 40 type 0

for {set id 1} {$id < $N} {incr id} {
	set dist 0.
	
#	while {$dist < 1. || pow($x-20., 2)+pow($y-20., 2)+pow($z-20., 2) < 11.*11.} {}
	while { $dist < 1 || [ constraint mindist_position $x $y $z ] < 1 } {
		set x [expr rand()*$box_l]
		set y [expr rand()*$box_l]
		set z [expr rand()*$box_l]
		
		part $id pos $x $y $z type 0
		set dist [analyze mindist]
	}
	
	puts "placed particle $id"
}

inter 0 0 lennard-jones 1. 1. 1.1225 0.25 0

prepare_vmd_connection "vmd" 2000

for {set i 0} {1} {incr i} {
	if {$i % 10 == 0} {
		lbfluid print vtk velocity velocity_$i.vtk	
		puts -nonewline "integrating $i"
		puts "th step"
		flush stdout
	}
	
  integrate 10
  imd positions
}
