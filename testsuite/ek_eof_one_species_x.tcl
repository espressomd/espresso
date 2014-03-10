# Copyright (C) 2011,2012,2013 The ESPResSo project
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

source "tests_common.tcl"

require_feature "ELECTROKINETICS"
require_feature "EK_BOUNDARIES"

set nn [format %02d [setmd n_nodes]]

puts "###############################################################"
puts "#         Testcase ek_eof_one_species_x.tcl running on        #"
puts "#                           $nn nodes                          #"
puts "###############################################################\n"

################################################################################
#                              Set up the System                               # 
################################################################################

# Set the slit pore geometry the width is the non-periodic part of the geometry
# the padding is used to ensure that there is no field inside outside the slit

set box_x 6
set box_y 6
set width 50

set padding 6
set box_z [expr $width+2*$padding]

setmd box_l $box_x $box_y $box_z

# Set the electrokinetic parameters

set agrid 1.0
set dt [expr 1.0/7.0]
set force 0.13
set sigma -0.05
set viscosity_dynamic 2.3
set friction 4.3
set temperature 2.9
set bjerrum_length 0.47

set temperature_LB [expr $agrid*$agrid/(3.0*$dt*$dt)]
set kB_LB 1.0
set cs_squared [expr (1.0/3.0)*($agrid*$agrid/($dt*$dt))]

# Set the simulation parameters

setmd time_step $dt
setmd skin 0.1
thermostat off
set integration_length 20000

# Set up the charged and neutral species

set density_neutral 33.4
set density_charged [expr -2.0*double($sigma)/double($width)]
set viscosity [expr $viscosity_dynamic*($density_charged + $density_neutral)]
set valency 1.0

# Set up the (LB) electrokinetics fluid

electrokinetics agrid $agrid lb_density [expr $density_neutral + $density_charged] viscosity $viscosity_dynamic friction $friction T $temperature bjerrum_length $bjerrum_length

electrokinetics 1 density $density_charged D 0.3 valency $valency ext_force $force 0 0
electrokinetics 2 density $density_neutral D 0.3 valency 0

# Set up the walls confining the fluid

electrokinetics boundary charge_density [expr $sigma/$agrid] rhomboid corner 0 0 [expr $padding-$agrid] c 0 0 $agrid b $box_x 0 0 a 0 $box_y 0 direction outside
electrokinetics boundary charge_density [expr $sigma/$agrid] rhomboid corner 0 0 [expr $padding+$width] c 0 0 $agrid b $box_x 0 0 a 0 $box_y 0 direction outside

# Set up the charged boundaries 

electrokinetics boundary charge_density 0.0 wall normal 0 0 1 d $padding 0 0 direction outside
electrokinetics boundary charge_density 0.0 wall normal 0 0 -1 d -[expr $padding+$width] 0 0 direction outside

# Integrate the system

integrate $integration_length

################################################################################
#                              Analyse the system                              # 
################################################################################

# Calculate the inverse length Xi, which is a combination of various
# constants (xi = q C / 2 kb T), with C a constant that needs to be
# solved for, or equivalently, xi needs to be solved for

# root finding function

proc solve {xi d bjerrum_length sigma valency } {
  set pi [expr {acos(-1.0)}]
  set el_char 1.0
  return [expr $xi*tan($xi*$d/2.0) + 2.0*$pi*$bjerrum_length*$sigma*$valency/$el_char ]
}

# initial parameters for bisection scheme

set pi [expr {acos(-1.0)}]
set size [expr $pi/(2.0*$width)]

set pnt0 0.0
set pntm [expr $pnt0 + $size]
set pnt1 [expr $pnt0 + 1.9*$size]

# the bisection scheme

set tol 1.0e-08
while { $size > $tol } {

  set val0 [solve $pnt0 $width $bjerrum_length $sigma $valency]
  set val1 [solve $pnt1 $width $bjerrum_length $sigma $valency]
  set valm [solve $pntm $width $bjerrum_length $sigma $valency]

  if { $val0 < 0.0 && $val1 > 0.0 } {
    if { $valm < 0.0 } {
      set pnt0 $pntm
      set size [expr $size/2.0]
      set pntm [expr $pnt0 + $size]
    } else {
      set pnt1 $pntm
      set size [expr $size/2.0]
      set pntm [expr $pnt1 - $size]
    }
  } elseif { $val0 > 0.0 && $val1 < 0.0 } {
    if { $valm < 0.0 } {
      set pnt1 $pntm
      set size [expr $size/2.0]
      set pntm [expr $pnt1 - $size]
    } else {
      set pnt0 $pntm
      set size [expr $size/2.0]
      set pntm [expr $pnt0 + $size]
    }
  } else {
    error_exit "Bisection method fails:\nTuning of domain boundaries may be required."
  }
}

# obtain the desired xi value

set xi $pntm

# function to calculate the density

proc density {x xi valency bjerrum_length} {
  set pi [expr {acos(-1.0)}]
  set kb 1.0
  return [expr ($xi*$xi)/(2.0*$pi*$bjerrum_length*$valency*$valency*cos($xi*$x)*cos($xi*$x)) ]
}

# function to calculate the velocity

proc velocity {x xi d valency bjerrum_length force viscosity} {
  set pi [expr {acos(-1.0)}]
  set el_char 1.0
  return [expr ($force)*log(cos($xi*$x)/cos($xi*$d/2.0))/(2.0*$pi*$bjerrum_length*$viscosity*$valency*$valency*$el_char*$el_char) ]
}

# function to calculate the xz component of the stress tensor

proc stress_tensor {x xi valency bjerrum_length force} {
  set pi [expr {acos(-1.0)}]
  set el_char 1.0
  return [expr $force*$xi*tan($xi*$x)/(2.0*$pi*$bjerrum_length*$valency*$valency*$el_char*$el_char) ]
}

# compare the various quantities to the analytic results

set total_velocity_difference 0.0
set total_density_difference 0.0
set total_stress_difference_xx 0.0
set total_stress_difference_yy 0.0
set total_stress_difference_zz 0.0
set total_stress_difference_xx_yy 0.0
set total_stress_difference_xx_zz 0.0
set total_stress_difference_yy_zz 0.0
set total_stress_difference_xy 0.0
set total_stress_difference_yz 0.0
set total_stress_difference_xz 0.0

for {set i 0} {$i < [expr $box_z/$agrid]} {incr i} {
  
  if {[expr $i*$agrid] >= $padding && [expr $i*$agrid] < [expr $box_z - $padding] } {
    set xvalue [expr $i*$agrid - $padding]
    set position [expr $i*$agrid - $padding - $width/2.0 + $agrid/2.0]

    # density
    set measured_density [electrokinetics 1 node [expr int($box_x/(2*$agrid))] [expr int($box_y/(2*$agrid))] $i print density]
    set calculated_density [density $position $xi $valency $bjerrum_length]
    set density_difference [expr 2.0*abs($measured_density - $calculated_density)/(abs($measured_density) + abs($calculated_density))]
    set total_density_difference [expr $total_density_difference + $density_difference ]

    # velocity
    set measured_velocity [lindex [lbnode [expr int($box_x/(2*$agrid))] [expr int($box_y/(2*$agrid))] $i print u] 0]
    set calculated_velocity [velocity $position $xi $width $valency $bjerrum_length $force $viscosity]
    set velocity_difference [expr 2.0*abs($measured_velocity - $calculated_velocity)/(abs($measured_velocity) + abs($calculated_velocity))]
    set total_velocity_difference [expr $total_velocity_difference + $velocity_difference ]

    # diagonal stress tensor

    set density_offset [electrokinetics 1 node [expr int($box_x/(2*$agrid))] [expr int($box_y/(2*$agrid))] [expr int($box_z/(2*$agrid))] print density]
    set measured_density [electrokinetics 1 node [expr int($box_x/(2*$agrid))] [expr int($box_y/(2*$agrid))] $i print density]
    set measured_density [expr $measured_density - $density_offset]

    set stress_0_offset [lindex [lbnode [expr int($box_x/(2*$agrid))] [expr int($box_y/(2*$agrid))] [expr int($box_z/(2*$agrid))] print pi_neq] 0]
    set measured_stress_0 [lindex [lbnode [expr int($box_x/(2*$agrid))] [expr int($box_y/(2*$agrid))] $i print pi_neq] 0]
    set measured_stress_0 [expr $measured_stress_0 - $stress_0_offset]

    set stress_2_offset [lindex [lbnode [expr int($box_x/(2*$agrid))] [expr int($box_y/(2*$agrid))] [expr int($box_z/(2*$agrid))] print pi_neq] 2]
    set measured_stress_2 [lindex [lbnode [expr int($box_x/(2*$agrid))] [expr int($box_y/(2*$agrid))] $i print pi_neq] 2]
    set measured_stress_2 [expr $measured_stress_2 - $stress_2_offset]

    set stress_5_offset [lindex [lbnode [expr int($box_x/(2*$agrid))] [expr int($box_y/(2*$agrid))] [expr int($box_z/(2*$agrid))] print pi_neq] 5]
    set measured_stress_5 [lindex [lbnode [expr int($box_x/(2*$agrid))] [expr int($box_y/(2*$agrid))] $i print pi_neq] 5]
    set measured_stress_5 [expr $measured_stress_5 - $stress_5_offset]    

    set stress_difference_xx [expr abs( $measured_stress_0 - $cs_squared*$density_offset*$agrid*$agrid )]
    set stress_difference_yy [expr abs( $measured_stress_2 - $cs_squared*$density_offset*$agrid*$agrid )]
    set stress_difference_zz [expr abs( $measured_stress_5 - $cs_squared*$density_offset*$agrid*$agrid )]

    set total_stress_difference_xx [expr $total_stress_difference_xx + $stress_difference_xx ]
    set total_stress_difference_yy [expr $total_stress_difference_yy + $stress_difference_yy ]
    set total_stress_difference_zz [expr $total_stress_difference_zz + $stress_difference_zz ]

    # diagonal stress tensor comparison
    set measured_stress_0 [lindex [lbnode [expr int($box_x/(2*$agrid))] [expr int($box_y/(2*$agrid))] $i print pi_neq] 0]
    set measured_stress_2 [lindex [lbnode [expr int($box_x/(2*$agrid))] [expr int($box_y/(2*$agrid))] $i print pi_neq] 2]
    set measured_stress_5 [lindex [lbnode [expr int($box_x/(2*$agrid))] [expr int($box_y/(2*$agrid))] $i print pi_neq] 5]
    set stress_difference_xx_yy [expr 2.0*abs($measured_stress_0 - $measured_stress_2)/(abs($measured_stress_0) + abs($measured_stress_2))]
    set stress_difference_yy_zz [expr 2.0*abs($measured_stress_5 - $measured_stress_2)/(abs($measured_stress_5) + abs($measured_stress_2))]
    set stress_difference_xx_zz [expr 2.0*abs($measured_stress_0 - $measured_stress_5)/(abs($measured_stress_0) + abs($measured_stress_5))]
    set total_stress_difference_xx_yy [expr $total_stress_difference_xx_yy + $stress_difference_xx_yy ]
    set total_stress_difference_yy_zz [expr $total_stress_difference_yy_zz + $stress_difference_yy_zz ]
    set total_stress_difference_xx_zz [expr $total_stress_difference_xx_zz + $stress_difference_xx_zz ]

    # xy component stress tensor
    set measured_stress [lindex [lbnode [expr int($box_x/(2*$agrid))] [expr int($box_y/(2*$agrid))] $i print pi_neq] 1]
    set calculated_stress 0.0
    set stress_difference_xy [expr abs($measured_stress - $calculated_stress)]
    set total_stress_difference_xy [expr $total_stress_difference_xy + $stress_difference_xy ]

    # yz component stress tensor
    set measured_stress [lindex [lbnode [expr int($box_x/(2*$agrid))] [expr int($box_y/(2*$agrid))] $i print pi_neq] 4]
    set calculated_stress 0.0
    set stress_difference_yz [expr abs($measured_stress - $calculated_stress)]
    set total_stress_difference_yz [expr $total_stress_difference_yz + $stress_difference_yz ]

    # xz component stress tensor
    set measured_stress [lindex [lbnode [expr int($box_x/(2*$agrid))] [expr int($box_y/(2*$agrid))] $i print pi_neq] 3]
    set calculated_stress [stress_tensor $position $xi $valency $bjerrum_length $force]
    set stress_difference_xz [expr 2.0*abs($measured_stress - $calculated_stress)/(abs($measured_stress) + abs($calculated_stress))]
    set total_stress_difference_xz [expr $total_stress_difference_xz + $stress_difference_xz ]
  }
}

set total_density_difference [expr $agrid*$total_density_difference/$box_z]
set total_velocity_difference [expr $agrid*$total_velocity_difference/$box_z]
set total_stress_difference_xx [expr $agrid*$total_stress_difference_xx/$box_z]
set total_stress_difference_yy [expr $agrid*$total_stress_difference_yy/$box_z]
set total_stress_difference_zz [expr $agrid*$total_stress_difference_zz/$box_z]
set total_stress_difference_xx_yy [expr $agrid*$total_stress_difference_xx_yy/$box_z]
set total_stress_difference_yy_zz [expr $agrid*$total_stress_difference_yy_zz/$box_z]
set total_stress_difference_xx_zz [expr $agrid*$total_stress_difference_xx_zz/$box_z]
set total_stress_difference_xy [expr $agrid*$total_stress_difference_xy/$box_z]
set total_stress_difference_yz [expr $agrid*$total_stress_difference_yz/$box_z]
set total_stress_difference_xz [expr $agrid*$total_stress_difference_xz/$box_z]

puts "Density deviation: $total_density_difference"
puts "Velocity deviation: $total_velocity_difference\n"
puts "Stress deviation xx component: $total_stress_difference_xx"
puts "Stress deviation yy component: $total_stress_difference_yy"
puts "Stress deviation zz component: $total_stress_difference_zz"
puts "Stress deviation xy component: $total_stress_difference_xy"
puts "Stress deviation yz component: $total_stress_difference_yz"
puts "Stress deviation xz component: $total_stress_difference_xz\n"
puts "Stress deviation between xx and yy: $total_stress_difference_xx_yy"
puts "Stress deviation between yy and zz: $total_stress_difference_yy_zz"
puts "Stress deviation between xx and zz: $total_stress_difference_xx_zz\n"
puts "NB. The stress on the diagonal of the tensor is only isotropic"
puts "    in equilibrium. However, it is not isotropic for the LB."
puts "    The anisotropic part relaxes towards isotropic, but it"
puts "    is not instantaneously isotropic. The elements on the"
puts "    diagonal must therefore be different.\n"

if { $total_density_difference > 2.5e-03 } {
  error_exit "Density accuracy not achieved"
}
if { $total_velocity_difference > 5.0e-03 } {
  error_exit "Velocity accuracy not achieved"
}
if { $total_stress_difference_xx > 1.0e-02 } {
  error_exit "Difference xx to yy component too large"
}
if { $total_stress_difference_yy > 1.0e-02 } {
  error_exit "Difference yy to zz component too large"
}
if { $total_stress_difference_zz > 1.0e-02 } {
  error_exit "Difference xx to zz component too large"
}
if { $total_stress_difference_xy > 5.0e-06 } {
  error_exit "Pressure accuracy xy component not achieved"
}
if { $total_stress_difference_yz > 5.0e-06 } {
  error_exit "Pressure accuracy yz component not achieved"
}
if { $total_stress_difference_xz > 5.0e-03 } {
  error_exit "Pressure accuracy xz component not achieved"
}
if { $total_stress_difference_xx_yy > 7.5e-03 } {
  error_exit "Difference xx to yy component too large"
}
if { $total_stress_difference_xx_yy < 2.5e-03 } {
  error_exit "Difference xx to yy component too small"
}
if { $total_stress_difference_yy_zz > 5.0e-06 } {
  error_exit "Difference yy to zz component too large"
}
if { $total_stress_difference_xx_zz > 7.5e-03 } {
  error_exit "Difference xx to zz component too large"
}
if { $total_stress_difference_xx_zz < 2.5e-03 } {
  error_exit "Difference xx to zz component too small"
}

exit 0
