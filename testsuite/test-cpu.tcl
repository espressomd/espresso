# the filenames must be set

set vmd "n"

setmd time_step 0.1;    # timestep for the integrator
setmd skin 0.2;
thermostat off;
set pipeX 100;
set pipeY 20;
set pipeZ 20;

setmd box_l $pipeX $pipeY $pipeZ;  # rectangular channel is defined

set originX 10;       # the initial coordinates of the center
set originY 10;       # of the immersed object
set originZ 10; 

set stretchX 1.0;     # the immersed object will be scaled 
set stretchY 1.0;     # by these factors in correspomding directions
set stretchZ 1.0;

set tmp 3
	# I will create a cuboidal mesh with tmpxtmpxtmp nodes.
for {set i 0} { $i < $tmp} {incr i} {
	for {set j 0} { $j < $tmp} {incr j} {
		for {set k 0} { $k < $tmp} {incr k} {
			part [expr $tmp*$tmp*$i + $tmp*$j + $k] pos [expr $originX + $i] [expr $originY + $j] [expr $originZ + $k] type 0;
		}
	}
}

inter 0 fene 1.0 3.0 1.0;
	# creating horizontal bonds
for {set i 0} { $i < [expr $tmp - 1]} {incr i} {
	for {set j 0} { $j < $tmp} {incr j} {
		for {set k 0} { $k < $tmp} {incr k} {
			part [expr $tmp*$tmp*$i + $tmp*$j + $k] bond 0 [expr $tmp*$tmp*($i + 1) + $tmp*$j + $k]
		}
	}
}	

	# creating vertical bonds
for {set i 0} { $i < $tmp} {incr i} {
	for {set j 0} { $j < [expr $tmp - 1]} {incr j} {
		for {set k 0} { $k < $tmp} {incr k} {
			part [expr $tmp*$tmp*$i + $tmp*$j + $k] bond 0 [expr $tmp*$tmp*$i + $tmp*($j + 1) + $k]
		}
	}
}	

	# creating front-back bonds
for {set i 0} { $i < $tmp} {incr i} {
	for {set j 0} { $j < $tmp} {incr j} {
		for {set k 0} { $k < [expr $tmp - 1]} {incr k} {
			part [expr $tmp*$tmp*$i + $tmp*$j + $k] bond 0 [expr $tmp*$tmp*$i + $tmp*$j + $k + 1]
		}
	}
}	



cellsystem domain_decomposition -no_verlet_list;      
lbfluid grid 1 dens 1.0 visc 1.5 tau 0.1 friction 0.5;
                           

if { $vmd == "y" } {
    prepare_vmd_connection simEspresso 3000 1
                     #visualization
    exec sleep 2   
    imd positions
}


# main iteration loop

set cycle 0 
while { $cycle<200} {
	puts "$cycle";
    if { $vmd == "y"} { imd positions};

	for {  set i 0 } { $i < [expr $tmp*$tmp*$tmp] } { incr i } {
			puts "[part $i print pos]"
		}


  # setting the constant velocity
  # of the fluid on the left side of the md_box
  for { set i 0 } { $i < 1} { incr i } {
    for { set j 0 } { $j < 20 } { incr j } {
      for { set k 0 } { $k < 20 } { incr k } {
        lbnode $i $j $k set u 0.5 0.0 0.0;
      }
    }
  }
  integrate 1;
  incr cycle;
}
