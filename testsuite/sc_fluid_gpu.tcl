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
set box_lz  1.0
set tstep 1.0
setmd skin 0.2
setmd time_step $tstep
setmd box_l  $box_lx $box_ly $box_lz
setmd periodic 1 1 0

# SC fluid parameters
set coupl 1.5  
set max   2.3 
set min   0.06
set avg   1.0
set visc 0.1
set mob 1e-1

if { [ catch {

lbfluid gpu dens $avg $avg mobility $mob  visc $visc $visc  friction 1.0 1.0   agrid 1.0 tau 1.0  sc_coupling 0.0 $coupl 0.0 
thermostat off 

for { set x  0  } { $x < $box_lx  } { incr x } { 
    for {set y 0} { $y < $box_ly } {incr y }  {
      for {set z 0} { $z < $box_lz } {incr z }  { 
	lbnode $x $y $z set rho [expr $min+(($max-$min)*(tanh(($box_lx/4-sqrt( ($x-$box_lx/2.)*($x-$box_lx/2.)  ))/2.)+1))/2. ]  [expr $max+(($min-$max)*(tanh(($box_lx/4-sqrt( ($x-$box_lx/2.)*($x-$box_lx/2.)   ))/2.)+1))/2. ] 
      }
    }
}

integrate 1000
set v {0 0 0}
set vnorm {0 0 0}
for {set z 0} { $z < $box_lz } {incr z }  { 
  for {set y 0} { $y < $box_ly } {incr y }  {
      for { set x  0  } { $x < $box_lx  } { incr x } { 
           set J [lbnode $x $y $z print v]
	   set v [vecadd $v $J]
	   set sqv ""
           foreach c $J { lappend sqv [expr $c*$c] }	   
	   set vnorm [vecadd $vnorm $sqv]
      }
    }
}
set average_momentum_per_node [vecscale [expr 1./(32*32)] $v]
set total_absolute_momentum_per_node  ""
foreach c [vecscale [expr 1./(32*32)]  $vnorm  ] { lappend total_absolute_momentum_per_node [expr sqrt($c)]}
#puts "$average_momentum_per_node  $total_absolute_momentum_per_node"
foreach c $average_momentum_per_node {
      if { [expr abs($c) > 1e-5 ]  } { 
            error "average momentum per node too large ([expr abs($c)])"
      }
}
foreach c $total_absolute_momentum_per_node {
      if { [expr abs($c) < 1e-20 ]  || [expr abs($c) > 1e-4 ]  } { 
           lbfluid print vtk density "relaxA.vtk" "relaxB.vtk"
           lbfluid print vtk velocity "relaxV.vtk" 
      }
      if { [expr abs($c) < 1e-20 ]  } { 
	   error "average momentum per node too close to zero, this could mean that the coupling is off, check the relax*.vtk files"
      }
      if { [expr abs($c) > 1e-4 ]  } { 
            error "average total absolute momentum per node too large ([expr abs($c)]), check the relax*.vtk files"
      }
}

} res ] } { 
    error_exit $res
}

exit 0

