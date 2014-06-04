# Copyright (C) 2012,2013 The ESPResSo project
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

#    TEST DESCRIPTION
#
#    This test case loads the model of a blood cell, 
#    puts it into a fluid, sets the nonzero fluid velocity 
#    on the left inflow and lets the blood cell flows for 100 timesteps.
#    
#    In the beginning, the configuration is loaded 
#    from object_in_fluid_system.data.init
#    The model consists of a triangular mesh with 400 nodes.
#
#    After 100 timesteps, the positions, velocities and forces of the 400 particles 
#    are stored in arrays POS, VEL and FOR.
# 
#    Then, the reference configuration is loaded 
#    from object_in_fluid_system.data.final
#    and a check is performed whether computed configuration
#    stored in FOR, VEL and POS corresponds to the reference configuration. 

source "tests_common.tcl"

require_feature "AREA_FORCE_GLOBAL"
require_feature "VOLUME_FORCE"
require_feature "LB_GPU"
require_feature "SHANCHEN" off

require_max_nodes_per_side 2

puts "------------------------------------------------"
puts "- Testcase object_in_fluid_gpu.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "------------------------------------------------"

set vmd "n"

set tcl_precision 15
set tolerance 1e-4

setmd time_step 0.1
setmd skin 0.2
thermostat off

proc read_data {file} {
    set f [open $file "r"]
    while {![eof $f]} { blockfile $f read auto}
    close $f
}

proc write_data_init {file} {
    set f [open $file "w"]
    blockfile $f write variable box_l
    blockfile $f write particles {id type mol pos}
    blockfile $f write interactions
    blockfile $f write bonds all
    close $f
}

proc write_data_final {file} {
    set f [open $file "w"]
    blockfile $f write particles {id type mol pos v f}
    close $f
}

if { [catch {
	
	read_data "object_in_fluid_system-init.data"

	lbfluid gpu grid 1 dens 1.0 visc 1.5 tau 0.1 friction 0.5
		                           
	if { $vmd == "y" } {
	    prepare_vmd_connection simEspresso 3000 1 
	    exec sleep 2   
	    imd positions
	}
	
	
	# main iteration loop
	
	set cycle 0 
	while { $cycle < 20 } {
	    puts -nonewline "time step $cycle/20\r"; flush stdout
	    if { $vmd == "y"} { imd positions}
	
	
	    # set the constant velocity
	    # of the fluid on the left side of the md_box
	    for { set i 0 } { $i < 1} { incr i } {
		for { set j 0 } { $j < 20 } { incr j } {
		    for { set k 0 } { $k < 20 } { incr k } {
			lbnode $i $j $k set u 0.5 0.0 0.0
		    }
		}
	    }

	    integrate 1

	    incr cycle
	}
	
	set generate_new_data 0
	# Here, you write new reference configuration in case you uncomment the next line
  # set generate_new_data 1
	if { $generate_new_data == 1} {
	    write_data_final "object_in_fluid_system-final.data"
	    exit
	}
	
	# store computed values for velocities, positions and forces
	for { set i 0 } { $i <= [setmd max_part] } { incr i } {
	    set POS($i) [part $i pr pos]
	    set VEL($i) [part $i pr v]
	    set FOR($i) [part $i pr f]
	}
	
	# load reference configuration
	read_data "object_in_fluid_system-final.data"
	
	set diffPOS 0.0
	set diffVEL 0.0
	set diffFOR 0.0
	
	for { set i 0 } { $i <= [setmd max_part] } { incr i } {
	    set tmp [part $i pr pos]
	    set Ax [lindex $tmp 0]
	    set Ay [lindex $tmp 1]
	    set Az [lindex $tmp 2]
	    # stores the vector of the reference position for $i-th particle
	    set Bx [lindex $POS($i) 0]
	    set By [lindex $POS($i) 1]
	    set Bz [lindex $POS($i) 2]
	    # stores the vector of the computed position for $i-th particle
	    set diffPOS [expr $diffPOS + sqrt(($Ax-$Bx)*($Ax-$Bx) + ($Ay-$By)*($Ay-$By) + ($Az-$Bz)*($Az-$Bz))]

	    set tmp [part $i pr v]
	    set Ax [lindex $tmp 0]
	    set Ay [lindex $tmp 1]
	    set Az [lindex $tmp 2]
	    # stores the vector of the reference velocity for $i-th particle
	    set Bx [lindex $VEL($i) 0]
	    set By [lindex $VEL($i) 1]
	    set Bz [lindex $VEL($i) 2]
	    # stores the vector of the computed velocity for $i-th particle
	    set diffVEL [expr $diffVEL + sqrt(($Ax-$Bx)*($Ax-$Bx) + ($Ay-$By)*($Ay-$By) + ($Az-$Bz)*($Az-$Bz))]

	    set tmp [part $i pr f]
	    set Ax [lindex $tmp 0]
	    set Ay [lindex $tmp 1]
	    set Az [lindex $tmp 2]
	    # stores the vector of the reference force for $i-th particle
	    set Bx [lindex $FOR($i) 0]
	    set By [lindex $FOR($i) 1]
	    set Bz [lindex $FOR($i) 2]
	    # stores the vector of the computed force for $i-th particle
	    set diffFOR [expr $diffFOR + sqrt(($Ax-$Bx)*($Ax-$Bx) + ($Ay-$By)*($Ay-$By) + ($Az-$Bz)*($Az-$Bz))]
	}

	puts "difference between the reference configuration and the computed configuration: "
	puts "		positions: $diffPOS"
	puts "		velocities: $diffVEL"
	puts "		forces: $diffFOR"
	if { $diffPOS > $tolerance || $diffVEL > $tolerance || $diffFOR > $tolerance } {
	    error "A difference occured between the reference configuration and the computed configuration"
	}
    } res ] } {
    error_exit $res
}

exit 0
