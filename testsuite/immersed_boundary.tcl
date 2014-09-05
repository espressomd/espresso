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
#    puts it into a fluid, sets the shear or poiseuille flow
#    on the left inflow and lets the blood cell flows for 100 timesteps.
#    
#    In the beginning, the configuration is loaded 
#    from immersed_boundary_system.data.init
#    The model consists of a triangular mesh with 642 nodes.
#
#    After 100 timesteps, the positions, velocities and forces of the 642 particles 
#    are stored in arrays POS, VEL and FOR.
# 
#    Then, the reference configuration is loaded 
#    from immersed_boundary_system.data.final
#    and a check is performed whether computed configuration
#    stored in FOR, VEL and POS corresponds to the reference configuration. 

source "tests_common.tcl"
source "immersed_boundary_system-tools.tcl"

require_feature "IMMERSED_BOUNDARY"
require_feature "VIRTUAL_SITES_IMMERSED_BOUNDARY"
require_feature "STRETCHING_FORCE_IMMERSED_BOUNDARY"
require_feature "VOLUME_CONSERVATION_IMMERSED_BOUNDARY"
require_feature "LB"
require_feature "LB_BOUNDARIES"

require_max_nodes_per_side 40

puts "------------------------------------------------"
puts "- Testcase immersed_boundary.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "------------------------------------------------"

##Test Parameters##

set tcl_precision 15
set tolerance 1e-8

proc read_data {file} {
    set f [open $file "r"]
    while {![eof $f]} { blockfile $f read auto}
    close $f
}

proc write_data_init {file} {
    set f [open $file "w"]
    blockfile $f write variable box_l
    blockfile $f write particles {id pos v f}
    close $f
}

proc write_data_final {file} {
    set f [open $file "w"]
    blockfile $f write particles {id pos v f}
    close $f
}


# Box geometry
set boxx  20
set boxy  20
set boxz  20

# Fix physical quantities (SI)

# kinematic viscosity 
# -> if visc is given here time step has to be set later
set nu_SI [expr 1e-6] 
# density of fluid
set rho_SI 1e3
# shear modulus
set ks_SI 5e-6
# bending
set kb_SI 1e-19
# grid size
set gridsize_SI 1e-7
# thermal energy
set kbT_SI [expr 1.38e-23 * 0]


# Particle parameters:

# radius:
# hydrodynamic radius, used to calculate friction constant. Diameter should be around grid constant
set rad_SI [expr 0.5*$gridsize_SI] 
# friction constant: is only set here for consitancy
set zeta_SI [expr 6*3.142*$nu_SI*$rho_SI*$rad_SI]
# # drag force on particle 
set part_force_SI [expr 0.01*$zeta_SI]
# particle volume 
set vol_part_SI [expr 4/3*3.142*$rad_SI*$rad_SI*$rad_SI]
# mass 
# set mass_SI [expr $rho_SI*$vol_part_SI]

########### Four possible unit systems ###################################

# # Unit system I

# # a priori fixed scales
#set rho0 $rho_SI
#set nu0 $nu_SI
#set L0 $deltaX_SI

# # derived scales
#set E0 [expr $rho0*$nu0*$nu0*$L0]
#set ks0 [expr $E0/($L0*$L0)]
#set V0 [expr $nu0/$L0]
#set kb0 [expr $E0]

# # Unit system II

# # a priori fixed scales
#set rho0 $rho_SI
#set V0 $V_SI
#set L0 $deltaX_SI

# # derived scales
#set nu0 [expr $L0*$V0]
#set E0 [expr $rho0*$nu0*$nu0*$L0]
#set ks0 [expr $E0/($L0*$L0)]
#set kb0 [expr $E0]

# Unit system  III

# fixed scales
set rho0 1e3
set L0 1e-7
# be careful energy scale is not set to 0°c!
set E0 [expr 1.38e-23 * 273]

# derived scales
set nu0 [expr sqrt($E0/($rho0*$L0))]
# set V0 $nu0/$L0
set ks0 [expr $E0/($L0*$L0) ]
set kb0 $E0
set m0 [expr $rho0*$L0*$L0*$L0]


# # Unit system  IV
# 
# # fixed scales  
# set rho0 $rho_SI
# set L0 $deltaX_SI
# set ks0 $ks_SI
# 
# # derived scales
# set E0 [expr $ks0*$L0*$L0]
# set nu0 [expr sqrt($E0/($rho0*$L0))]
# set V0 $nu0/$L0 
# set kb0 $E0

############### Calculate quantities in Espresso units ###############################
set pi 3.1415926535897931
set G 0.05
set r 4.9856099295288183473
set H 2
set rho [expr $rho_SI / $rho0]
set nu [expr $nu_SI / $nu0]
set kbT [expr $kbT_SI / $E0]
set kS [expr $ks_SI / $ks0]
set kB [expr $kb_SI / $kb0]
set kC $kS

set shear_rate [expr ($G * $kS) / ($nu * $r)]
set u [expr $shear_rate * ($boxz - (($boxz - $H)/2))]

# Hagen–Poiseuille equation #
set LP [expr $boxz - $H]
set QP [expr $shear_rate*($boxz - $H) * $pi * $LP*$LP]
set PP [expr (8 * $nu * $boxz * $QP) / ($pi * $LP*$LP*$LP*$LP)]
set FP [expr $PP / $boxz]

# set V [expr $V_SI /  $V0]
set gridsize [expr $gridsize_SI / $L0]

# particle scales
#radius
set rad [expr $rad_SI / $L0]
# volume
set vol_part_SI [expr 4/3*3.142*$rad*$rad*$rad]
# mass 
set mass 1
# friction constant: is calculated new from viscosity and paricle radius
set zeta [expr 6*3.142*$nu*$rho*$rad]
# drag force
set part_force [expr 0.01*$zeta]

########################################################################
# Numerical parameters

# maximum strech length for bonds
set maxStretch 3
# timestep calculated from condition that relaxation parameter tau in LB has to be 1 (or similar)
set taurelax 1
set dt_lb [expr ($taurelax - 0.5)*$gridsize*$gridsize/(3*$nu)]
set dt_md [expr $dt_lb]

# Dimensionless numbers following Le
# both in espresso and SI units for consistency check
# particle radius needs to be compatible with node file
# lengthscale of simulation
set L [expr $boxz-1] 
set L_SI [expr $L * $L0]
# dynamic viscosity
set mu [expr $nu * $rho]
set mu_SI [expr $nu_SI * $rho_SI]
# thermal mean squared velocity
set v_mean [expr sqrt(3 * $kbT / $mass)]
set mass_SI [expr $mass*$m0]
set v_mean_SI [expr sqrt(3 * $kbT_SI / $mass_SI)]

# Renolds number calculated with meansquared velocity
set Re [expr $rad*$v_mean / $nu]
set Re_SI [expr $rad_SI*$v_mean_SI / $nu_SI]

# "stiffness ratio" - probably should be about one??
# set Sr [expr $kC/$kS]
# set Sr_SI [expr $kc_SI / $ks_SI]

# output
puts "nu = $nu"
puts "rho = $rho"
puts "Re = $Re     Re_SI = $Re_SI"
puts "dt_lb = $dt_lb"
puts "dt_md = $dt_md"
puts "kS = $kS"
puts "kB = $kB"
puts "kC = $kC"
# puts "kbT = $kbT0"
puts "zeta = $zeta"
puts "Radius = $r"
puts "Shear Rate =  [expr ($G * $kS)] / [expr ($nu * $r)] * $H =  [expr ( ($G * $kS) / ($nu * $r) )]"
puts "u = [expr $shear_rate * ($boxz - ($H/2))]"


#########################################################################
# SETUP 
#########################################################################

# setting Boxlength
setmd box_l $boxx $boxy $boxz

# setting integration parameters
# skin for verlet list
setmd skin 0.1
# timestep
setmd time_step $dt_md
# coose how good no slip condition is: 0 non, 1 better
# setmd sequ 1

setmd warnings 0

#setting themostat
thermostat lb $kbT

############## INITIALIZE FLOW  ###############

set flow 0

if {$flow == 0} {

set flowType "Shear flow with volume conservation"
# setting up the fluid with or without using gpu
lbfluid agrid $gridsize dens $rho visc $nu tau $dt_lb friction $zeta ext_force 0 0 0
# walls located at z=1 and z=boxz-1 in x-y-plane
lbboundary wall dist [expr $H/2] normal 0. 0. 1. velocity $u 0 0 
lbboundary wall dist [expr -$boxz + ($H/2)] normal 0. 0. -1. velocity -$u 0 0

} elseif {$flow == 1} { 

set flowType "Hagen-Pouiessele Flow with volume conservation"
# setting up the fluid with or without using gpu
lbfluid agrid $gridsize dens $rho visc $nu tau $dt_lb friction $zeta ext_force $FP 0 0
# walls located at z=1 and z=boxz-1 in x-y-plane
lbboundary wall dist [expr $H/2] normal 0. 0. 1. velocity 0 0 0 
lbboundary wall dist [expr -$boxz + ($H/2)] normal 0. 0. -1. velocity 0 0 0

}

################## ADD CELLS #############
set interCount 0
set cellCount 0
set nodeCount 0
set numNodesPerCell 0

set numType 0

# number of cells
set numCells 1
# setmd vescnum $numCells

addCell [expr $boxx*0.5] [expr $boxy*0.5] [expr $boxz*0.5] $numType
# addCell 10 10 10 $numType
incr numType


if { [catch {

    set write_init_data 0
    # Here, you write the initial reference configuration in case you uncomment the next line
    # set write_init_data 1
    if { $write_init_data == 1} {
	write_data_init "immersed_boundary_system-init.data"
	exit
    }

	
	# main iteration loop
	
	set cycle 0 
	while { $cycle < 100 } {
	    puts -nonewline "time step $cycle/100\r"; flush stdout
       
	    integrate 1

	    incr cycle
	}
	
	set generate_new_data 0
	# Here, you write new reference configuration in case you uncomment the next line
  # set generate_new_data 1
	if { $generate_new_data == 1} {
	    write_data_final "immersed_boundary_system-generate.data"
	    exit
	}
	
	# store computed values for velocities, positions and forces
	for { set i 0 } { $i <= [setmd max_part] } { incr i } {
	    set POS($i) [part $i pr pos]
	    set VEL($i) [part $i pr v]
	    set FOR($i) [part $i pr f]
	}
	
	# load reference configuration
	read_data "immersed_boundary_system-final.data"
	
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
	    error "A difference occured between the reference configuration and the computed configuration that is higher than the allowed tolerance value"
	}
    } res ] } {
    error_exit $res
}

puts "OK: particle properties correspondence within a $flowType after $cycle steps"


exit 0
