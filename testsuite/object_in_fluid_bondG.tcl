# Copyright (C) 2010,2011,2012,2013 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#   Max-Planck-Institute for Polymer Research, Theory Group
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
require_feature "LB"


puts "------------------------------------------------"
puts "- Testcase object_in_fluid.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "------------------------------------------------"

set createPart TMP/createPart;

set bondS TMP/bondsStretching
set bondB TMP/bondsBending
set bondAlocal TMP/bondsAreaLocal
set bondAglobal TMP/bondsAreaGlobal
set bondV TMP/bondsVolume
set bondVA TMP/bondsVolumeAreaGlobal
set partS TMP/partStretching
set partB TMP/partBending
set partAlocal TMP/partAreaLocal
set partAglobal TMP/partAreaGlobal
set partV TMP/partVolume
set partVA TMP/partVolumeAreaGlobal


set vmd "n"

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

	# 	Here you can write new initial configuration - by changing generate_new_data from 0 to 1
	# 	For this, you need to supply the meshfiles: one for nodes and one for triangles
	
	set generate_new_data 0	
	if { $generate_new_data == 1} {
		puts "generating new data"
		set fileNodes object_in_fluid_system.data.nodes
		set fileTriangles object_in_fluid_system.data.triangles
		setmd box_l 100 20 20

		set ks 0.05;          # stretching of the cell
		set kb 0.1;          # bending
		set kal 0.1;         # local area preservation
		set kag 0.01;         # global area preservation
		set kv 10.;           # volume preservation
		
		set massPart 1.0;     # setting the mass of praticles
		
		set nnode 400;          # number of IB points on the surface
		                      # of the immersed object
		                      # for example with an red blood cell
		                      # set 400
		set ntriangle 796;      # number of triangles at the surface
		                      # of the immersed object
		                      # for example with an red blood cell
		                      # set 796
		                      
		set originX 10.1;       # the initial coordinates of the center
		set originY 10.1;       # of the immersed object
		set originZ 10.1; 
		
		set stretchX 1.0;     # the immersed object will be scaled 
		set stretchY 1.0;     # by these factors in correspomding directions
		set stretchZ 1.0;
		
		set rotateX 0.;       # Rotation by specified angles around 
		set rotateY 0.;       # X axis, Y axis and Z axis. Angles given
		set rotateZ 0.;       # in radians. rotateX=Pi/2 rotates the object 
		                      # by 90 degrees with the axes of rotation x 
		                      # such that vector 0,1,0 changes to 0,0,1 and
		                      # 0,0,1 changes to 0,-1,0
		                      
		set typePart 0;       # each immersed object must have 
		set molPart 0;        # different type and mol ID
			
		set firstBondId 0;    # when adding more immersed objects, 
		set firstPartId 0;    # these parameters need to be set properly
		
		
		# the following script generates the source files
		# for particles, interactions and bonds
		
		exec ./bondGenerator $fileNodes $fileTriangles $nnode \
		$ntriangle $rotateX $rotateY $rotateZ $originX $originY $originZ \
		$stretchX $stretchY $stretchZ $createPart$molPart \
		$typePart $molPart $firstPartId  $bondS$molPart \
		$bondB$molPart $bondAlocal$molPart $bondAglobal$molPart \
		$bondV$molPart $bondVA$molPart $partS$molPart $partB$molPart \
		$partAlocal$molPart $partAglobal$molPart $partV$molPart $partVA$molPart\
		$ks $kb $kal $kag $kv $firstBondId $massPart;
		
		# create particles and bonds from source files
		
		source $createPart$molPart; 
		source $bondB$molPart; 
		source $partB$molPart; 

		write_data_init "object_in_fluid_bondG_system.data.init"
	} else {
		read_data "object_in_fluid_bondG_system.data.init"
		invalidate_system
	}
	
	cellsystem domain_decomposition -no_verlet_list 
	lbfluid grid 1 dens 1.0 visc 1.5 tau 0.1 friction 0.5	
	
	                           
	if { $vmd == "y" } {
	    prepare_vmd_connection simEspresso 3000 1 
	    exec sleep 2   
	    imd positions
	}
	
	
	# main iteration loop
	
	set cycle 0 
	while { $cycle<100 } {
		puts "$cycle";
	    if { $vmd == "y"} { imd positions};
	
	
	  # set the constant velocity
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
	
	# Here, you write new reference configuration in case you have chosen to generate new data
	#
	if { $generate_new_data == 1} {
		write_data_final "object_in_fluid_bondG_system.data.final"
		exit
	}
	
	# store computed values for velocities, positions and forces
	for { set i 0 } { $i <= [setmd max_part] } { incr i } {
		set POS($i) [part $i pr pos]
		set VEL($i) [part $i pr v]
		set FOR($i) [part $i pr f]
	}
	
	# load reference configuration
	 read_data "object_in_fluid_bondG_system.data.final"
	
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
	if { $diffPOS > 0.0 || $diffPOS > 0.0 || $diffPOS > 0.0 } {
		puts "A difference occured between the reference configuration and the computed configuration: "
		puts "		positions: $diffPOS"
		puts "		velocities: $diffVEL"
		puts "		forces: $diffFOR"
		error "A difference occured between the reference configuration and the computed configuration"
    }
} res ] } {
    error_exit $res
}

exit 0
