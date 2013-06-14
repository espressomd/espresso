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

require_feature "LB_GPU"
require_feature "LB_BOUNDARIES_GPU"
require_feature "EXTERNAL_FORCES"

puts "---------------------------------------------------------------"
puts "- Testcase lb_planar_embedded_particles_gpu.tcl running on [format %02d [setmd n_nodes]] nodes"
puts "---------------------------------------------------------------"

# Here we test different features of the LB subsytem in a planar slit geometry.
# The test is meant to check the agreement of the LB calculation with
# analytically known results under different parameter combinations.

# The box length is set to 12, because this allows many convenient lattice constants
set l 12
setmd box_l $l $l $l
setmd time_step 0.01

set agrid 0.75
set visc 7.
set rho 2.
set tau 0.04

setmd skin [ expr 0.4*$agrid ]

lbfluid gpu agrid $agrid visc $visc dens $rho tau $tau
lbfluid friction 1
thermostat lb 0

for {set i 0} {$i < 10} {incr i} {
    part $i pos $i 0 0
}

# Wall positions are set to $agrid and $l-$agrid, leaving one layer of boundary nodes 
# on each side of the box
#
# We start with couette flow setting the velocity on boundary 2 to $v_couette

set v_couette 0.000001

set wall_normal_direction 1
set normal1 [ list 0 0 0 ]
lset normal1 $wall_normal_direction 1
set normal2 [ list 0 0 0 ]
lset normal2 $wall_normal_direction -1

set couette_flow_direction 2
set v_boundary [ list 0 0 0 ]
lset v_boundary $couette_flow_direction $v_couette

lbboundary wall normal [ lindex $normal1 0 ]  [ lindex $normal1 1 ]  [ lindex $normal1 2 ] dist $agrid
lbboundary wall normal [ lindex $normal2 0 ]  [ lindex $normal2 1 ]  [ lindex $normal2 2 ] dist [ expr -$l +$agrid ] \
  velocity [ lindex $v_boundary 0 ] [ lindex $v_boundary 1 ] [ lindex $v_boundary 2 ] 

set dist [ expr $l - 2*$agrid ]

integrate 2000

set accuracy_u 0.
set meanabs_u 0.
set accuracy_p 0.
set meanabs_p 0.

for { set i 2 } { $i < int(floor($l/$agrid))-2 } { incr i } {
  set pos [ expr ($i+0.5)*$agrid ]

  set u_theory [ expr $v_couette*($pos-$agrid)/$dist ]
  set p_xy_theory [ expr -$v_couette/$dist*$rho*$visc ]
  ## We go though the box in normal direction
  set m [ expr $i * [ lindex $normal1 0 ] ]
  set n [ expr $i * [ lindex $normal1 1 ] ]
  set o [ expr $i * [ lindex $normal1 2 ] ]
  set u_measured [ lindex [ lbnode $m $n $o print u ] $couette_flow_direction ]
  if { ($couette_flow_direction == 0 && $wall_normal_direction == 1) || ($couette_flow_direction == 1 && $wall_normal_direction == 0) } {
    set p_xy_measured  [ lindex [ lbnode $m $n $o print pi_neq ] 1 ]
  }
  if { ($couette_flow_direction == 0 && $wall_normal_direction == 2) || ($couette_flow_direction == 2 && $wall_normal_direction == 0) } {
    set p_xy_measured  [ lindex [ lbnode $m $n $o print pi_neq ] 3 ]
  }
  if { ($couette_flow_direction == 1 && $wall_normal_direction == 2) || ($couette_flow_direction == 2 && $wall_normal_direction == 1) } {
    set p_xy_measured  [ lindex [ lbnode $m $n $o print pi_neq ] 4 ]
  }
  # invert sign of p_xy
#  puts "$p_xy_measured"
  set accuracy_u [ expr $accuracy_u + abs($u_theory - $u_measured) ]
  set meanabs_u [ expr $meanabs_u + abs($u_theory) ]
  set accuracy_p [ expr $accuracy_u + abs($p_xy_theory - $p_xy_measured ) ]
  set meanabs_p [ expr $meanabs_p + abs($p_xy_theory) ]
}

set couette_u_accuracy [ expr $accuracy_u / $meanabs_u ]
set couette_p_accuracy  [ expr $accuracy_p / $meanabs_p ]

puts "Couette flow result:"
puts "flow accuracy $couette_u_accuracy"
puts "pressure accuracy $couette_p_accuracy"
puts "----------"

# Now we add a force density in normal direction, and compress the flow.
# We expect that the pressure gradient exactly cancels the force density:
# grad p = -f_hydrostatic
set f_hydrostatic 0.1

set fx [ expr $f_hydrostatic*[ lindex $normal1 0 ] ]
set fy [ expr $f_hydrostatic*[ lindex $normal1 1 ] ]
set fz [ expr $f_hydrostatic*[ lindex $normal1 2 ] ]

lbfluid gpu agrid $agrid visc $visc dens $rho friction 1 tau $tau ext_force $fx $fy $fz
integrate 1000

# Equation of state: p = rho*c_2**2
set p_center [ expr $rho * 1./3. * $agrid *$agrid/$tau/$tau ]

set accuracy_u 0.
set meanabs_u 0.
set accuracy_p 0.
set meanabs_p 0.
for { set i 2 } { $i < int(floor($l/$agrid))-2 } { incr i } {
  set pos [ expr ($i+0.5)*$agrid ]
  set p_neq_xx_theory [ expr $f_hydrostatic*($pos - $l/2) ]
  set p_xx_theory [ expr $p_neq_xx_theory + $p_center ]

  set m [ expr $i * [ lindex $normal1 0 ] ]
  set n [ expr $i * [ lindex $normal1 1 ] ]
  set o [ expr $i * [ lindex $normal1 2 ] ]

  if { $wall_normal_direction == 0 } {
    set p_xx_measured  [ lindex [ lbnode $m $n $o print pi ] 0 ]
    set p_neq_xx_measured  [ lindex [ lbnode $m $n $o print pi_neq ] 0 ]
  }
  if { $wall_normal_direction == 1 } {
    set p_xx_measured  [ lindex [ lbnode $m $n $o print pi ] 2 ]
    set p_neq_xx_measured  [ lindex [ lbnode $m $n $o print pi_neq ] 2 ]
  }
  if { $wall_normal_direction == 2 } {
    set p_xx_measured  [ lindex [ lbnode $m $n $o print pi ] 5 ]
    set p_neq_xx_measured  [ lindex [ lbnode $m $n $o print pi_neq ] 5 ]
  }
  set accuracy_p [ expr $accuracy_u + abs($p_neq_xx_theory - $p_neq_xx_measured ) ]
  set meanabs_p [ expr $meanabs_p + abs($p_neq_xx_theory) ]
}
set hydrostatic_p_accuracy  [ expr $accuracy_p / $meanabs_p ]
puts "Hydrostatic test result:"
puts "pressure accuracy $hydrostatic_p_accuracy"
puts "-------------"

# Now we add a force density in the direction of the Couette flow
# We'd like to switch of the Couette flow, but changing the velocity of the BCs is not implemented
# yet. So we just make a much larger force density.
# 
# We expect a flow profile u=f/eta/2*(x-agrid)*($l-$agrid-$x)
# and the pressure tensor has to be the derivative of that: p_xy = f*(x-$l/2)

set f_body 0.1

set f_body_vec [ list 0 0 0 ]
lset f_body_vec $couette_flow_direction $f_body

lbfluid gpu agrid $agrid visc $visc dens $rho friction 1 tau $tau ext_force \
  [ lindex $f_body_vec 0 ] [ lindex $f_body_vec 1 ] [ lindex $f_body_vec 2 ]

integrate 2000

set accuracy_u 0.
set meanabs_u 0.
set accuracy_p 0.
set meanabs_p 0.
for { set i 2 } { $i < int(floor($l/$agrid))-2 } { incr i } {
  set pos [ expr ($i+0.5)*$agrid ]

  set u_theory [ expr  $f_body/$visc/$rho/2*($pos-$agrid)*($l-$agrid-$pos) ]
  set p_xy_theory [ expr $f_body * ($pos-$l/2) ]
  ## We go though the box in normal direction
  set m [ expr $i * [ lindex $normal1 0 ] ]
  set n [ expr $i * [ lindex $normal1 1 ] ]
  set o [ expr $i * [ lindex $normal1 2 ] ]
  set u_measured [ lindex [ lbnode $m $n $o print u ] $couette_flow_direction ]
  if { ($couette_flow_direction == 0 && $wall_normal_direction == 1) || ($couette_flow_direction == 1 && $wall_normal_direction == 0) } {
    set p_xy_measured  [ lindex [ lbnode $m $n $o print pi_neq ] 1 ]
  }
  if { ($couette_flow_direction == 0 && $wall_normal_direction == 2) || ($couette_flow_direction == 2 && $wall_normal_direction == 0) } {
    set p_xy_measured  [ lindex [ lbnode $m $n $o print pi_neq ] 3 ]
  }
  if { ($couette_flow_direction == 1 && $wall_normal_direction == 2) || ($couette_flow_direction == 2 && $wall_normal_direction == 1) } {
    set p_xy_measured  [ lindex [ lbnode $m $n $o print pi_neq ] 4 ]
  }
#  puts "$pos $u_theory $u_measured $p_xy_theory $p_xy_measured"
  set accuracy_u [ expr $accuracy_u + abs($u_theory - $u_measured) ]
  set meanabs_u [ expr $meanabs_u + abs($u_theory) ]
  set accuracy_p [ expr $accuracy_u + abs($p_xy_theory - $p_xy_measured ) ]
  set meanabs_p [ expr $meanabs_p + abs($p_xy_theory) ]
}
set poisseuille_u_accuracy [ expr $accuracy_u / $meanabs_u ]
set poisseuille_p_accuracy [ expr $accuracy_p / $meanabs_p ]

puts "Poisseuille flow result:"
puts "flow accuracy $poisseuille_u_accuracy"
puts "pressure accuracy $poisseuille_p_accuracy"
puts "----------"

if { $couette_u_accuracy > 1e-5 } {
  error_exit "Couette flow accuracy not achieved"
}
if { $couette_p_accuracy > 1e-5 } {
  error_exit "Couette pressure accuracy not achieved"
}
if { $hydrostatic_p_accuracy > 1e-3 } {
  error_exit "hydrostatic pressure accuracy not achieved"
}
if { $poisseuille_u_accuracy > 3e-2 } {
  error_exit "Poisseuille flow accuracy not achieved"
}
if { $poisseuille_p_accuracy > 1e-2 } {
  error_exit "Poisseuille pressure accuracy not achieved"
}

exit 0
