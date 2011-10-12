set box_l 40
set temp 1
set N 1000


setmd skin 0.1
setmd box_l $box_l $box_l $box_l
setmd time_step 0.01

thermostat langevin $temp 0.1

constraint rhomboid corner 20 20 20 a 5 0 0 b 0 5 0 c 0 0 5 direction outside type 0
#constraint sphere center 20 20 20 radius 10 direction outside type 0

part 0 pos 40 40 40 type 0

for {set id 1} {$id < $N} {incr id} {
	set dist 0.
	
	while {$dist < 1. || pow($x-20., 2)+pow($y-20., 2)+pow($z-20., 2) < 11.*11.} {
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
	if {$i % 100 == 0} {
		puts -nonewline "integrating $i"
		puts "th step"
		flush stdout
	}
	
  integrate 1
  imd positions
}
