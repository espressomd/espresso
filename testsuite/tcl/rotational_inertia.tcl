# Copyright (C) 2011,2012,2013,2014,2015,2016 The ESPResSo project
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
source "tests_common.tcl"

require_feature "MASS"
require_feature "ROTATIONAL_INERTIA"

proc define_rotation_matrix { part } {
    
    # 3*3 matrix
    set A {list 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0}
    
    for {set k 0} {$k<4} {incr k} {
        set quat($k) [lindex [part $part print quat] $k]
    }
    
    for {set k 0} {$k<3} {incr k} {
        set quatu($k) [lindex [part $part print quatu] $k]
    }
    
    set q0q0 [expr $quat(0) * $quat(0)]
    set q1q1 [expr $quat(1) * $quat(1)]
    set q2q2 [expr $quat(2) * $quat(2)]
    set q3q3 [expr $quat(3) * $quat(3)]
        
    
    lset A [expr 0 + 3*0] [expr $q0q0 + $q1q1 - $q2q2 - $q3q3]
    lset A [expr 1 + 3*1] [expr $q0q0 - $q1q1 + $q2q2 - $q3q3]
    lset A [expr 2 + 3*2] [expr $q0q0 - $q1q1 - $q2q2 + $q3q3]
    
    lset A [expr 0 + 3*1] [expr 2*($quat(1)*$quat(2) + $quat(0)*$quat(3))]
    lset A [expr 0 + 3*2] [expr 2*($quat(1)*$quat(3) - $quat(0)*$quat(2))]
    lset A [expr 1 + 3*0] [expr 2*($quat(1)*$quat(2) - $quat(0)*$quat(3))]
    
    lset A [expr 1 + 3*2] [expr 2*($quat(2)*$quat(3) + $quat(0)*$quat(1))]
    lset A [expr 2 + 3*0] [expr 2*($quat(1)*$quat(3) + $quat(0)*$quat(2))]
    lset A [expr 2 + 3*1] [expr 2*($quat(2)*$quat(3) - $quat(0)*$quat(1))]
    
    return $A
}

proc convert_vec_body_to_space { part vec } {
    set A [define_rotation_matrix $part]
    set vec_lab {list 0.0 0.0 0.0}

    lset vec_lab 0 [expr [expr [lindex $A [expr 0 + 3*0]]*[lindex $vec 0]] + [expr [lindex $A [expr 1 + 3*0]]*[lindex $vec 1]] + [expr [lindex $A [expr 2 + 3*0]]*[lindex $vec 2]]]
    lset vec_lab 1 [expr [expr [lindex $A [expr 0 + 3*1]]*[lindex $vec 0]] + [expr [lindex $A [expr 1 + 3*1]]*[lindex $vec 1]] + [expr [lindex $A [expr 2 + 3*1]]*[lindex $vec 2]]]
    lset vec_lab 2 [expr [expr [lindex $A [expr 0 + 3*2]]*[lindex $vec 0]] + [expr [lindex $A [expr 1 + 3*2]]*[lindex $vec 1]] + [expr [lindex $A [expr 2 + 3*2]]*[lindex $vec 2]]]
    
   
    return $vec_lab
}

# Angular momentum
proc L_body { part } {
    set L_body_return [list 0.0 0.0 0.0]
    for {set k 0} {$k<3} {incr k} {
        lset L_body_return $k [expr [lindex [part $part print omega_body] $k] * [lindex [part $part print rinertia] $k]]
    }
    
    return $L_body_return
}

setmd skin 0
part 0 pos 0 0 0

#Inertial motion around the stable and unstable axes

thermostat off
set T "0.0 0.0 0.0"
# Anisotropic inertial moment. Stable axes correspond to J[1] and J[2].
# The unstable axis corresponds to J[0]. These values relation is J[1] < J[0] < J[2].
set J "5 0.5 18.5"

# Validation of J[1] stability
setmd time_step 0.0006
# Stable omega component should be larger than other components.
set stable_omega 57.65
part 0 omega_body 0.15 $stable_omega -0.043 ext_torque [lindex $T 0] [lindex $T 1] [lindex $T 2] rinertia [lindex $J 0] [lindex $J 1] [lindex $J 2]

# Angular momentum
set L_0_body [list 0.0 0.0 0.0]
set L_0_body [L_body 0]
set L_0_lab [convert_vec_body_to_space 0 $L_0_body]

for {set i 0} {$i <100} {incr i} {
   set L_body [L_body 0]
   set L_lab [convert_vec_body_to_space 0 $L_body]
   for {set k 0} {$k<3} {incr k} {
       if { abs([lindex $L_lab $k] - [lindex $L_0_lab $k]) > 4E-3 } {
        error_exit "Inertial motion around stable axis J1: Deviation in angular momentum is too large. Step $i, coordinate $k, expected [lindex $L_0_lab $k], got [lindex $L_lab $k]"
    }
   }
   set expected [expr $stable_omega]
   if { abs([lindex [part 0 print omega_body] 1] - $expected) > 4E-3 } {
      error_exit "Inertial motion around stable axis J1: Deviation in omega is too large. Step $i, coordinate 1, expected $expected, got [lindex [part 0 print omega_lab] 1]"
    }
  integrate 10
}

# Validation of J[2] stability
setmd time_step 0.01
set stable_omega 3.2
part 0 omega_body 0.011 -0.043 $stable_omega ext_torque [lindex $T 0] [lindex $T 1] [lindex $T 2] rinertia [lindex $J 0] [lindex $J 1] [lindex $J 2]

set L_0_body [L_body 0]
set L_0_lab [convert_vec_body_to_space 0 $L_0_body]

for {set i 0} {$i <100} {incr i} {
   set L_body [L_body 0]
   set L_lab [convert_vec_body_to_space 0 $L_body]
   for {set k 0} {$k<3} {incr k} {
       if { abs([lindex $L_lab $k] - [lindex $L_0_lab $k]) > 4E-3 } {
        error_exit "Inertial motion around stable axis J2: Deviation in angular momentum is too large. Step $i, coordinate $k, expected [lindex $L_0_lab $k], got [lindex $L_lab $k]"
    }
   }
   set expected [expr $stable_omega]
   if { abs([lindex [part 0 print omega_body] 2] - $expected) > 4E-3 } {
      error_exit "Inertial motion around stable axis J2: Deviation in omega too large. Step $i, coordinate 2, expected $expected, got [lindex [part 0 print omega_lab] 2]"
    }
  integrate 10
}

# Validation of J[0]
setmd time_step 0.001
# Unstable omega component should be larger than other components.
set unstable_omega 5.76
part 0 omega_body $unstable_omega -0.043 0.15 ext_torque [lindex $T 0] [lindex $T 1] [lindex $T 2] rinertia [lindex $J 0] [lindex $J 1] [lindex $J 2]

set L_0_body [L_body 0]
set L_0_lab [convert_vec_body_to_space 0 $L_0_body]

for {set i 0} {$i <100} {incr i} {
   set L_body [L_body 0]
   set L_lab [convert_vec_body_to_space 0 $L_body]
   for {set k 0} {$k<3} {incr k} {
       if { abs([lindex $L_lab $k] - [lindex $L_0_lab $k]) > 4E-3 } {
        error_exit "Inertial motion around unstable axis J0: Deviation in angular momentum is too large. Step $i, coordinate $k, expected [lindex $L_0_lab $k], got [lindex $L_lab $k]"
    }
   }
  integrate 10
}

part delete

exit 0
