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

require_feature "LB"
require_feature "EXTERNAL_FORCES"

puts "---------------------------------------------------------------"
puts "- Testcase lb_fluid_coupling.tcl running on [format %02d [setmd n_nodes]] nodes"
puts "---------------------------------------------------------------"

set tcl_precision 10

# box related properties
set length 12.0
setmd box_l $length $length $length
set box_vol [expr $length*$length*$length]
setmd periodic 1 1 1

# set the time step and skin
set tstep 0.1
setmd time_step $tstep
setmd skin 0.3

# set the external force in an irrational direction
set dragx 0.3308117918333205
set dragy 0.28623688095048233
set dragz 0.8992396823804915

set fdragx [expr -$dragx/$box_vol]
set fdragy [expr -$dragy/$box_vol]
set fdragz [expr -$dragz/$box_vol]

# set the lbfluid and thermostat
lbfluid cpu agrid 1 dens 1.0 visc 3.0 tau $tstep ext_force $fdragx $fdragy $fdragz friction 10.0
thermostat lb 0.0

# set the particle
part 0 pos [expr 0.5*$length] [expr 0.5*$length] [expr 0.5*$length] v 0.0 0.0 0.0 f 0.0 0.0 0.0 ext_force $dragx $dragy $dragz 

# get over the initial acceleration
integrate 200

# average terminal velocity to remove node dependence
set vsum 0.0
set count 0
for { set i 0 } { $i < 100 } { incr i } {

  integrate 5

  set pv [part 0 print v]
  set vel [expr sqrt([lindex $pv 0]*[lindex $pv 0] + [lindex $pv 1]*[lindex $pv 1] + [lindex $pv 2]*[lindex $pv 2]) ]
  set vsum [expr $vsum + $vel]
  incr count
}

# check for the right terminal velocity
set vel_works 0.1100128183

set difference [expr abs(($vsum/$count - $vel_works)/$vel_works)]

if { $difference > 1e-3 } {
  error_exit "Particle terminal velocity is wrong: coupling might be broken."
}

exit 0
