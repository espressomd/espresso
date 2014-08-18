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

require_feature "IMMERSED_BOUNDARY"
require_feature "VIRTUAL_SITES_IMMERSED_BOUNDARY"
require_feature "BENDING_FORCE_IMMERSED_BOUNDARY"
require_feature "EXTERNAL_FORCES"
require_feature "LB"
require_feature "LB_BOUNDARIES"

set error_margin 1e-5
set disp_margin 1e-1

set pi [expr atan(1) * 4]
proc vcrossp {a b} { 
    set v0 [expr ([lindex $a 1] * [lindex $b 2]) - ([lindex $a 2] * [lindex $b 1])]
    set v1 [expr ([lindex $a 2] * [lindex $b 0]) - ([lindex $a 0] * [lindex $b 2])]
    set v2 [expr ([lindex $a 0] * [lindex $b 1]) - ([lindex $a 1] * [lindex $b 0])]
    
    return [list $v0 $v1 $v2]

}

proc vdotp {a b} { 
    
    set sc [expr ([lindex $a 0] * [lindex $b 0]) + ([lindex $a 1] * [lindex $b 1]) + ([lindex $a 2] * [lindex $b 2])]
	    return $sc
}

proc vnor {a} {

    set mag [expr sqrt(([lindex $a 0] * [lindex $a 0]) + ([lindex $a 1] * [lindex $a 1]) + ([lindex $a 2] * [lindex $a 2]))]
   
    set nor0 [expr [lindex $a 0] / $mag]
    set nor1 [expr [lindex $a 1] / $mag]
    set nor2 [expr [lindex $a 2] / $mag]
    
    return [list $nor0 $nor1 $nor2]
    
}

puts "---------------------------------------------------------------"
puts "- Testcase triel.tcl running on  [format %02d [setmd n_nodes]] nodes: -"
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

# Verify that the triangle tries to maintain it's shape during $numSteps steps after a single x + 1.0 stretch to a vertex

# Set maximum bond stretch length
set maxStretch 2

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

set numSteps 100

part 0 pos 2 2 0 virtual 1
part 1 pos 4 2 0 virtual 1
part 2 pos 2 5 0 virtual 1
part 3 pos 0 2 0 virtual 1

inter 2 tribend 0 1 2 3 0 $kB $maxStretch
part 0 bond 2 1 2 3

puts [part 0 print pos ]
puts [part 1 print pos ]
puts [part 2 print pos ]
puts [part 3 print pos ]

set part0PosRef [part 0 print pos]
set part1PosRef [part 1 print pos]
set part2PosRef [part 2 print pos]
set part3PosRef [part 3 print pos]

part 1 pos [lindex $part1PosRef 0] [lindex $part1PosRef 1] [expr [lindex $part1PosRef 2] + 1.5]
part 3 pos [lindex $part3PosRef 0] [lindex $part3PosRef 1] [expr [lindex $part3PosRef 2] + 1.5]

integrate 0
if { ([part 0 print pos] != "2.0 2.0 0.0") || ([part 1 print pos] != "4.0 2.0 1.5")  || ([part 2 print pos] != "2.0 5.0 0.0") || ([part 3 print pos] != "0.0 2.0 1.5")} {
 error_exit  "Error: lbtracers incorrectly positioned: [part 0 print pos] [part 1 print pos] [part 2 print pos] [part 3 print pos]"
} else {
 puts "OK: Position of lbtracers"
}

set dx1 [vecsub $part0PosRef $part2PosRef]
set dx2 [vecsub $part1PosRef $part2PosRef]
set dx3 [vecsub $part3PosRef $part2PosRef]

set n1 [vcrossp $dx1 $dx2]
set n2 [vcrossp $dx1 $dx3]

lset n2 0 [expr [lindex $n2 0] * -1]
lset n2 1 [expr [lindex $n2 1] * -1]
lset n2 2 [expr [lindex $n2 2] * -1]

set n1 [vnor $n1]
set n2 [vnor $n2]

set sc [vdotp $n1 $n2]

set thetaRef [expr acos($sc)] 
set direc [vcrossp $n1 $n2]
set desc [vdotp $dx1 $direc]
if {$desc < 0} { 
    set thetaRef [expr (2 * $pi) - $thetaRef]
}
 
set error 0
puts "Checking "
for {set step 0} {$step < $numSteps} {incr step} {
    integrate 1

}

set part0PosComp [part 0 print pos]
set part1PosComp [part 1 print pos]
set part2PosComp [part 2 print pos]
set part3PosComp [part 3 print pos]

set dx1 [vecsub $part0PosComp $part2PosComp]
set dx2 [vecsub $part1PosComp $part2PosComp]
set dx3 [vecsub $part3PosComp $part2PosComp]

set n1 [vcrossp $dx1 $dx2]
set n2 [vcrossp $dx1 $dx3]

lset n2 0 [expr [lindex $n2 0] * -1]
lset n2 1 [expr [lindex $n2 1] * -1]
lset n2 2 [expr [lindex $n2 2] * -1]

set n1 [vnor $n1]
set n2 [vnor $n2]

set sc [vdotp $n1 $n2]
	
set thetaComp [expr acos($sc)] 
set direc [vcrossp $n1 $n2]
set desc [vdotp $dx1 $direc]

if {$desc < 0} { 
    set thetaComp [expr (2 * $pi) - $thetaComp]
}

if { $thetaRef  >= $thetaComp } {
    
    
    error_exit "Error: Bending angle between the two triangles not reduced. thetaRef: $thetaRef >= thetaComp: $thetaComp"
}

puts "OK: The two triangles try to maintain the reference bending angle during $numSteps steps after a single x + 1.5  stretch to two vertices in the z direction (Flat)"

part delete

# Verify that the triangle tries to return to it's bent shape

# Set maximum bond stretch length
set maxStretch 2

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

set numSteps 100

part 0 pos 2 2 0 virtual 1
part 1 pos 4 2 -1.5 virtual 1
part 2 pos 2 5 0 virtual 1
part 3 pos 0 2 -1.5 virtual 1

inter 2 tribend 0 1 2 3 0 $kB $maxStretch
part 0 bond 2 1 2 3

puts [part 0 print pos ]
puts [part 1 print pos ]
puts [part 2 print pos ]
puts [part 3 print pos ]

set part0PosRef [part 0 print pos]
set part1PosRef [part 1 print pos]
set part2PosRef [part 2 print pos]
set part3PosRef [part 3 print pos]

part 1 pos [lindex $part1PosRef 0] [lindex $part1PosRef 1] [expr [lindex $part1PosRef 2] + 1.5]
part 3 pos [lindex $part3PosRef 0] [lindex $part3PosRef 1] [expr [lindex $part3PosRef 2] + 1.5]

integrate 0
if { ([part 0 print pos] != "2.0 2.0 0.0") || ([part 1 print pos] != "4.0 2.0 0.0")  || ([part 2 print pos] != "2.0 5.0 0.0") || ([part 3 print pos] != "0.0 2.0 0.0")} {
 error_exit  "Error: lbtracers incorrectly positioned: [part 0 print pos] [part 1 print pos] [part 2 print pos] [part 3 print pos]"
} else {
 puts "OK: Position of lbtracers"
}

set dx1 [vecsub $part0PosRef $part2PosRef]
set dx2 [vecsub $part1PosRef $part2PosRef]
set dx3 [vecsub $part3PosRef $part2PosRef]

set n1 [vcrossp $dx1 $dx2]
set n2 [vcrossp $dx1 $dx3]

lset n2 0 [expr [lindex $n2 0] * -1]
lset n2 1 [expr [lindex $n2 1] * -1]
lset n2 2 [expr [lindex $n2 2] * -1]

set n1 [vnor $n1]
set n2 [vnor $n2]

set sc [vdotp $n1 $n2]

set thetaRef [expr acos($sc)] 
set direc [vcrossp $n1 $n2]
set desc [vdotp $dx1 $direc]
if {$desc < 0} { 
    set thetaRef [expr (2 * $pi) - $thetaRef]
}
 
set error 0
puts "Checking "
for {set step 0} {$step < $numSteps} {incr step} {
    integrate 1

}

set part0PosComp [part 0 print pos]
set part1PosComp [part 1 print pos]
set part2PosComp [part 2 print pos]
set part3PosComp [part 3 print pos]

set dx1 [vecsub $part0PosComp $part2PosComp]
set dx2 [vecsub $part1PosComp $part2PosComp]
set dx3 [vecsub $part3PosComp $part2PosComp]

set n1 [vcrossp $dx1 $dx2]
set n2 [vcrossp $dx1 $dx3]

lset n2 0 [expr [lindex $n2 0] * -1]
lset n2 1 [expr [lindex $n2 1] * -1]
lset n2 2 [expr [lindex $n2 2] * -1]

set n1 [vnor $n1]
set n2 [vnor $n2]

set sc [vdotp $n1 $n2]
	
set thetaComp [expr acos($sc)] 
set direc [vcrossp $n1 $n2]
set desc [vdotp $dx1 $direc]

if {$desc < 0} { 
    set thetaComp [expr (2 * $pi) - $thetaComp]
}

if { $thetaRef  <= $thetaComp } {
    
    
    error_exit "Error: Bending angle between the two triangles not increased. thetaRef: $thetaRef <= thetaComp: $thetaComp"
}

puts "OK: The two triangles try to maintain the reference bending angle during $numSteps steps after a single x + 1.5  stretch to two vertices in the z direction (Bent)"

exit 0
