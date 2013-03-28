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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 
#

### Measuring the force on a single sphere immersed in a fluid with
### fixed velocity boundary conditions created by two 
### walls at finite distance.
### The force is compared to th analytical result F=6 pi eta r v
### i.e. the stokes force on the particles.

source "tests_common.tcl"

require_feature "LB"
require_feature "LB_BOUNDARIES"

puts "---------------------------------------------------------------"
puts "- Testcase lb_stokes_sphere.tcl running on [format %02d [setmd n_nodes]] nodes"
puts "---------------------------------------------------------------"

set w 16
set l 16
setmd box_l [ expr $w+2 ] $l $l
setmd time_step 0.01
thermostat lb 0.
setmd skin 0.5
cellsystem domain_decomposition -no_verlet_list

set kinematic_visc 10.
set dens 1.
lbfluid cpu visc [ expr $kinematic_visc / $dens ] dens $dens friction 1. agrid 1.0 tau .01 

set v1 1.0
lbboundary wall normal -1. 0. 0. dist [ expr -(+0.5+$w) ] velocity 0.00 $v1 0.
lbboundary wall normal 1. 0. 0. dist 0.5 velocity 0. $v1 0.
set radius 2.5
lbboundary sphere center 8 8 8 radius $radius direction +1
integrate 1000
puts "[ lindex [ lbboundary force 0 ] 1 ] [ lindex [ lbboundary force 1 ] 1 ]"
set lbforce [ lindex [ lbboundary force 2 ] 1 ]
set twoRoverH [expr 2*$radius/$w]
set refforce [ expr 6*3.1415*$kinematic_visc*$radius*$v1 ]

#set refforce [ expr 6*3.1415*$kinematic_visc*$radius*$v1*(1-1.004*$twoRoverH+0.418*pow($twoRoverH,3)+0.210*pow($twoRoverH,4)-0.169*pow($twoRoverH,5)) ]
# Note the force comes for diffusion between plates from Balducci et al., Macromolecules 2006, 39, 6273-6281
set deviation [ expr abs( $lbforce -$refforce ) / $refforce ]
puts "The measured force is: $lbforce"
puts "The stokes force is: $refforce"
puts "deviation to stokes force is [ expr $deviation * 100 ] %"
if { $deviation > 0.2 } {
  error_exit "The deviation from Stokes law larger than expected"
}

exit 0
