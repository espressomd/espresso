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
#############################################################
#                                                           #
#  Test virutal sites                                       #
#                                                           #
#############################################################
source "tests_common.tcl"

require_feature "LBTRACERS"
require_feature "EXTERNAL_FORCES"
require_feature "LB"
require_feature "LB_BOUNDARY"

set error_margin 1e-5 

puts "---------------------------------------------------------------"
puts "- Testcase lbtracers.tcl running on 1 nodes"
puts "---------------------------------------------------------------"

# Box geometry
set boxx  10
set boxy  10
set boxz  10

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
set rho [expr $rho_SI / $rho0]
set nu [expr $nu_SI / $nu0]
set kbT [expr $kbT_SI / $E0]
set kS [expr $ks_SI / $ks0]
set kB [expr $kb_SI / $kb0]
set kC $kS

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

# timestep calculated from condition that relaxation parameter tau in LB has to be 1 (or similar)
set taurelax 1
set dt_lb [expr ($taurelax - 0.5)*$gridsize*$gridsize/(3*$nu)]
set dt_md [expr $dt_lb]
set numSteps 100

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

# output
puts "nu = $nu"
puts "rho = $rho"
puts "Re = $Re     Re_SI = $Re_SI"
puts "dt_lb = $dt_lb"
puts "dt_md = $dt_md"
puts "kS = $kS"
puts "kB = $kB"
puts "kC = $kC"
puts "kbT = $kbT"
puts "zeta = $zeta"

# setting Boxlength
setmd box_l $boxx $boxy $boxz

# setting integration parameters
# skin for verlet list
setmd skin 0.1
# timestep
setmd time_step $dt_md

setmd warnings 0

# setting up the fluid with or without using gpu
lbfluid agrid $gridsize dens $rho visc $nu tau $dt_lb friction $zeta ext_force 0 0 0
#setting themostat
thermostat lb $kbT

part 1 pos 5 5 6 virtual 
part 2 pos 5 5 4 virtual

puts [part 1 print pos ]
puts [part 2 print pos ]


integrate 0
if { ([part 1 print pos] != "5.0 5.0 6.0") || ([part 2 print pos] != "5.0 5.0 4.0")} {
 error_exit  "Error: lbtracers incorrectly positioned: [part 1 print pos] [part 2 print pos]"
} else {
 puts "OK: Position of lbtracers"
}

set error 0
puts "Checking "
for {set step 0} {$step < $numSteps} {incr step} {
integrate 1
# Verify that velocity of the lbtracers matches the velocity of fluid near their positions
set part1v [part 1 print v]
set part2v [part 2 print v]

set part1pos [part 1 print pos]
set part2pos [part 2 print pos]

set part1lbv [lbfluid print_interpolated_velocity [lindex $part1pos 0] [lindex $part1pos 1] [lindex $part1pos 2]]
set part2lbv [lbfluid print_interpolated_velocity [lindex $part2pos 0] [lindex $part2pos 2] [lindex $part2pos 2]]

if { [expr abs([lindex $part1v 0] - [lindex $part1lbv 0])] > $error_margin ||     
[expr abs([lindex $part1v 1] - [lindex $part1lbv 1])] > $error_margin || 
[expr abs([lindex $part1v 2] - [lindex $part1lbv 2])] > $error_margin} {

 error_exit "Error: Particle 1 velocity incorrect."
}

if { [expr abs([lindex $part2v 0] - [lindex $part2lbv 0])] > $error_margin ||     
[expr abs([lindex $part2v 1] - [lindex $part2lbv 1])] > $error_margin || 
[expr abs([lindex $part2v 2] - [lindex $part2lbv 2])] > $error_margin} {

 error_exit "Error: Particle 2 velocity incorrect."
}
}
puts "OK: Velocities of lbtracer particles matches the lb fluid velocity at near their positions"

# Testing wether periodic boundaries are handled correctly
part delete

# setting Boxlength
setmd box_l $boxx $boxy $boxz

# setting integration parameters
# skin for verlet list
setmd skin 0.1
# timestep
setmd time_step $dt_md

setmd warnings 0

# setting up the fluid with or without using gpu
lbfluid agrid $gridsize dens $rho visc $nu tau $dt_lb friction $zeta ext_force 1 0 0
#setting themostat
thermostat lb $kbT

part 1 pos 9 5 2 virtual 
part 2 pos 9 2 3 virtual

puts [part 1 print pos ]
puts [part 2 print pos ]


set error 0
puts "Checking "
for {set step 0} {$step < $numSteps} {incr step} {
integrate 1
# Verify that velocity of the lbtracers matches the velocity of fluid near their positions
set part1v [part 1 print v]
set part2v [part 2 print v]

set part1pos [part 1 print pos]
set part2pos [part 2 print pos]

set part1lbv [lbfluid print_interpolated_velocity [lindex $part1pos 0] [lindex $part1pos 1] [lindex $part1pos 2]]
set part2lbv [lbfluid print_interpolated_velocity [lindex $part2pos 0] [lindex $part2pos 2] [lindex $part2pos 2]]

if { [expr abs([lindex $part1v 0] - [lindex $part1lbv 0])] > $error_margin ||     
[expr abs([lindex $part1v 1] - [lindex $part1lbv 1])] > $error_margin || 
[expr abs([lindex $part1v 2] - [lindex $part1lbv 2])] > $error_margin} {

 error_exit "Error: Particle 1 velocity incorrect."
}

if { [expr abs([lindex $part2v 0] - [lindex $part2lbv 0])] > $error_margin ||     
[expr abs([lindex $part2v 1] - [lindex $part2lbv 1])] > $error_margin || 
[expr abs([lindex $part2v 2] - [lindex $part2lbv 2])] > $error_margin} {

 error_exit "Error: Particle 2 velocity incorrect."
}
}
puts "OK: Handling of periodic boundaries"


exit 0
