# Copyright (C) 2011,2012 The ESPResSo project
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
require_feature "EXTERNAL_FORCES"
require_feature "SHANCHEN"

puts "---------------------------------------------------------------"
puts "- Testcase sc_fluid_gpu.tcl running on [format %02d [setmd n_nodes]] nodes"
puts "---------------------------------------------------------------"

set tcl_precision 15
set box_lx  32.0
set box_ly  32.0
set box_lz  4.0
set tstep 1.0
setmd skin 0.0
setmd time_step $tstep
setmd box_l  $box_lx $box_ly $box_lz
setmd periodic 1 1 0

# SC fluid parameters
set coupl 1.5  
set max   5.3 
set min   0.06
set avg   1.0
set visc 0.1
set mob 1e-1

#part 0 pos [expr 3. * $box_lx / 2. - 0.75 ] [expr $box_ly / 2.] [expr $box_lz / 2.]  fix  type 0 
#part 1 pos [expr 3. * $box_lx / 2. + 0.75 ] [expr $box_ly / 2.] [expr $box_lz / 2.]  fix  type 0 
part 0 pos [expr 1.122462048 + 0.1 ] [expr $box_ly / 2.] [expr $box_lz / 2.]  fix  type 0 
part 1 pos [expr 0           ] [expr $box_ly / 2.] [expr $box_lz / 2.]  fix  type 0 
if { [ catch {
puts "initializing fluid"
lbfluid gpu dens $avg $avg mobility $mob  visc $visc $visc  friction 1.0 1.0   agrid 1.0 tau 100.0  sc_coupling 0.0 $coupl 0.0 
thermostat off 
puts "creating fluid"
for { set x  0  } { $x < $box_lx  } { incr x } { 
    for {set y 0} { $y < $box_ly } {incr y }  {
      for {set z 0} { $z < $box_lz } {incr z }  { 
	lbnode $x $y $z set rho [expr $min+(($max-$min)*(tanh(($box_lx/4-sqrt( ($x-$box_lx/2.)*($x-$box_lx/2.)  ))/2.)+1))/2. ]  [expr $max+(($min-$max)*(tanh(($box_lx/4-sqrt( ($x-$box_lx/2.)*($x-$box_lx/2.)   ))/2.)+1))/2. ] 
      }
    }
}
puts "placing particle"
inter 0 0 lennard-jones 1.0 1.0 2.0
inter 0 0 affinity  0.0 1.0
#puts "[inter]"
integrate 1 
puts "[part 0 print composition ] [part 1 print composition ] [part 0 print f]"
set force [lindex [part 0 print f] 0 ]
if {  $force > -0.00095 } {
	error_exit "force too large, affinity interaction might be broken"
}

} res ] } { 
    error_exit $res

}
exit 0

