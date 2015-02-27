# Copyright (C) 2011,2012,2013,2014 The ESPResSo project
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
puts "#   Testcase ek_eof_one_species_y_nonlinear.tcl running on    #"
puts "#                           $nn nodes                          #"
puts "###############################################################\n"

################################################################################
#                              Set up the System                               # 
################################################################################

# Set the slit pore geometry the width is the non-periodic part of the geometry
# the padding is used to ensure that there is no field inside outside the slit

set box_y 6
set box_z 6
set width 50

set padding 6
set box_x [expr $width+2*$padding]

setmd box_l $box_x $box_y $box_z

# Set the electrokinetic parameters

set agrid 0.5
set dt [expr 1.0/11.0]
set force 0.13
set sigma -0.03
set viscosity_kinematic 1.0
set friction 1.0
set temperature 2.3
set bjerrum_length 0.7

set temperature_LB [expr $agrid*$agrid/(3.0*$dt*$dt)]
set kB_LB 1.0
set cs_squared [expr (1.0/3.0)*($agrid*$agrid/($dt*$dt))]

# Set the simulation parameters

setmd time_step $dt
setmd skin 0.1
thermostat off
set integration_length 40000

# Set up the charged and neutral species

set density_water 26.15
set density_counterions [expr -2.0*double($sigma)/double($width)]
set valency 1.0

# Set up the (LB) electrokinetics fluid

electrokinetics agrid $agrid lb_density $density_water viscosity $viscosity_kinematic friction $friction T $temperature bjerrum_length $bjerrum_length stencil nonlinear

electrokinetics 1 density $density_counterions D 0.3 valency $valency ext_force 0 $force 0

# Set up the walls confining the fluid

electrokinetics boundary charge_density [expr $sigma/$agrid] rhomboid corner [expr $padding-$agrid] 0 0 b $agrid 0 0 a 0 $box_y 0 c 0 0 $box_z direction outside
electrokinetics boundary charge_density [expr $sigma/$agrid] rhomboid corner [expr $padding+$width] 0 0 b $agrid 0 0 a 0 $box_y 0 c 0 0 $box_z direction outside

# Set up the charged boundaries 

electrokinetics boundary charge_density 0.0 wall normal 1 0 0 d $padding 0 0 direction outside
electrokinetics boundary charge_density 0.0 wall normal -1 0 0 d -[expr $padding+$width] 0 0 direction outside

# Integrate the system

integrate $integration_length

################################################################################
#                              Analyze the system                              # 
################################################################################

# Calculate the inverse length xi, which is a combination of various
# constants (xi = zeC/2kBT), with C a constant that needs to be
# solved for, or equivalently, xi needs to be solved for

# root finding function

proc solve {xi d bjerrum_length sigma valency } {
  set pi [expr {acos(-1.0)}]
  set el_char 1.0
  return [expr $xi*tan($xi*$d/2.0) + 2.0*$pi*$bjerrum_length*$sigma/($valency*$el_char)]
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

proc density {x xi bjerrum_length} {
  set pi [expr {acos(-1.0)}]
  set kb 1.0
  return [expr ($xi*$xi)/(2.0*$pi*$bjerrum_length*cos($xi*$x)*cos($xi*$x)) ]
}

# function to calculate the velocity

proc velocity {x xi d bjerrum_length force viscosity_kinematic density_water} {
  set pi [expr {acos(-1.0)}]
  return [expr ($force)*log(cos($xi*$x)/cos($xi*$d/2.0))/(2.0*$pi*$bjerrum_length*$viscosity_kinematic*$density_water) ]
}

# function to calculate the xz component of the stress tensor

proc stress_tensor {x xi bjerrum_length force} {
  set pi [expr {acos(-1.0)}]
  return [expr $force*$xi*tan($xi*$x)/(4.0*$pi*$bjerrum_length) ]
}

# function to calculate the pressure

proc pressure {x xi bjerrum_length} {
  global box_x box_y box_z agrid
  set pi [expr {acos(-1)}]
  set offset [lindex [lbnode [expr int($box_x/(2*$agrid))] [expr int($box_y/(2*$agrid))] [expr int($box_z/(2*$agrid))] print pi_neq] 0]
  return [expr 3.0*$xi*$xi*tan($xi*$x)*tan($xi*$x)/(2.0*$pi*$bjerrum_length) + $offset]
}

# compare the various quantities to the analytic results

set total_velocity_difference 0.0
set total_density_difference 0.0
set total_stress_difference_xx 0.0
set total_stress_difference_yy 0.0
set total_stress_difference_zz 0.0
set total_stress_difference_xy 0.0
set total_stress_difference_yz 0.0
set total_stress_difference_xz 0.0

for {set i 0} {$i < [expr $box_x/$agrid]} {incr i} {

  if {[expr $i*$agrid] >= $padding && [expr $i*$agrid] < [expr $box_x - $padding] } {
    set xvalue [expr $i*$agrid - $padding]
    set position [expr $i*$agrid - $padding - $width/2.0 + $agrid/2.0]

    # density
    set measured_density [electrokinetics 1 node $i [expr int($box_y/(2*$agrid))] [expr int($box_z/(2*$agrid))] print density]
    set calculated_density [density $position $xi $bjerrum_length]
    set density_difference [expr abs($measured_density - $calculated_density)]
    set total_density_difference [expr $total_density_difference + $density_difference ]

    # velocity
    set measured_velocity [lindex [lbnode $i [expr int($box_y/(2*$agrid))] [expr int($box_z/(2*$agrid))] print u] 1]
    set calculated_velocity [velocity $position $xi $width $bjerrum_length $force $viscosity_kinematic $density_water]
    set velocity_difference [expr abs($measured_velocity - $calculated_velocity)]
    set total_velocity_difference [expr $total_velocity_difference + $velocity_difference ]

    # diagonal stress tensor
    set measured_stress_xx [lindex [lbnode $i [expr int($box_y/(2*$agrid))] [expr int($box_z/(2*$agrid))] print pi_neq] 0]
    set calculated_stress_xx [pressure $position $xi $bjerrum_length]
    set measured_stress_yy [lindex [lbnode $i [expr int($box_y/(2*$agrid))] [expr int($box_z/(2*$agrid))] print pi_neq] 2]
    set calculated_stress_yy [pressure $position $xi $bjerrum_length]
    set measured_stress_zz [lindex [lbnode $i [expr int($box_y/(2*$agrid))] [expr int($box_z/(2*$agrid))] print pi_neq] 5]
    set calculated_stress_zz [pressure $position $xi $bjerrum_length]

    set stress_difference_xx [expr abs($measured_stress_xx - $calculated_stress_xx)]
    set stress_difference_yy [expr abs($measured_stress_yy - $calculated_stress_yy)]
    set stress_difference_zz [expr abs($measured_stress_zz - $calculated_stress_zz)]

    set total_stress_difference_xx [expr $total_stress_difference_xx + $stress_difference_xx ]
    set total_stress_difference_yy [expr $total_stress_difference_yy + $stress_difference_yy ]
    set total_stress_difference_zz [expr $total_stress_difference_zz + $stress_difference_zz ]

    # xy component stress tensor
    set measured_stress_xy [lindex [lbnode $i [expr int($box_y/(2*$agrid))] [expr int($box_z/(2*$agrid))] print pi_neq] 1]
    set calculated_stress_xy [stress_tensor $position $xi $bjerrum_length $force]
    set stress_difference_xy [expr abs($measured_stress_xy - $calculated_stress_xy)]
    set total_stress_difference_xy [expr $total_stress_difference_xy + $stress_difference_xy ]

    # yz component stress tensor
    set measured_stress_yz [lindex [lbnode $i [expr int($box_y/(2*$agrid))] [expr int($box_z/(2*$agrid))] print pi_neq] 4]
    set calculated_stress_yz 0.0
    set stress_difference_yz [expr abs($measured_stress_yz - $calculated_stress_yz)]
    set total_stress_difference_yz [expr $total_stress_difference_yz + $stress_difference_yz ]

    # xz component stress tensor
    set measured_stress_xz [lindex [lbnode $i [expr int($box_y/(2*$agrid))] [expr int($box_z/(2*$agrid))] print pi_neq] 3]
    set calculated_stress_xz 0.0
    set stress_difference_xz [expr abs($measured_stress_xz - $calculated_stress_xz)]
    set total_stress_difference_xz [expr $total_stress_difference_xz + $stress_difference_xz ]
  }
}

set total_density_difference [expr $agrid*$total_density_difference/$box_x]
set total_velocity_difference [expr $agrid*$total_velocity_difference/$box_x]
set total_stress_difference_xx [expr $agrid*$total_stress_difference_xx/$box_x]
set total_stress_difference_yy [expr $agrid*$total_stress_difference_yy/$box_x]
set total_stress_difference_zz [expr $agrid*$total_stress_difference_zz/$box_x]
set total_stress_difference_xy [expr $agrid*$total_stress_difference_xy/$box_x]
set total_stress_difference_yz [expr $agrid*$total_stress_difference_yz/$box_x]
set total_stress_difference_xz [expr $agrid*$total_stress_difference_xz/$box_x]

puts "Density deviation: $total_density_difference"
puts "Velocity deviation: $total_velocity_difference\n"
puts "Stress deviation xx component: $total_stress_difference_xx"
puts "Stress deviation yy component: $total_stress_difference_yy"
puts "Stress deviation zz component: $total_stress_difference_zz"
puts "Stress deviation xy component: $total_stress_difference_xy"
puts "Stress deviation yz component: $total_stress_difference_yz"
puts "Stress deviation xz component: $total_stress_difference_xz\n"
puts "NB. The stress on the diagonal of the tensor is only isotropic"
puts "    in equilibrium. However, it is not isotropic for the LB."
puts "    The anisotropic part relaxes towards isotropic, but it"
puts "    is not instantaneously isotropic. The elements on the"
puts "    diagonal must therefore be different.\n"

if { $total_density_difference > 5.0e-07 } {
  error_exit "Density accuracy not achieved"
}
if { $total_velocity_difference > 5.0e-07 } {
  error_exit "Velocity accuracy not achieved"
}
#if { $total_stress_difference_xx > 5.0e-06 } { #TODO put this back in once the LB pressure output is fixed
#  error_exit "Difference xx component too large"
#}
#if { $total_stress_difference_yy > 7.5e-06 } {
#  error_exit "Difference yy component too large"
#}
#if { $total_stress_difference_zz > 5.0e-06 } {
#  error_exit "Difference zz component too large"
#}
if { $total_stress_difference_xy > 1.0e-06 } {
  error_exit "Stress accuracy xy component not achieved"
}
if { $total_stress_difference_yz > 5.0e-11 } {
  error_exit "Stress accuracy yz component not achieved"
}
if { $total_stress_difference_xz > 5.0e-11 } {
  error_exit "Stress accuracy xz component not achieved"
}

exit 0
