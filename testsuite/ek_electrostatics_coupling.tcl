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
require_feature "EK_ELECTROSTATIC_COUPLING"


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
set viscosity_kinematic 2.3
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

# Set up the charged and neutral species

set density_water 26.15
set density_counterions [expr -2.0*double($sigma)/double($width)]
set valency 1.0

# Set up the (LB) electrokinetics fluid

electrokinetics agrid $agrid lb_density $density_water viscosity $viscosity_kinematic friction $friction T $temperature bjerrum_length $bjerrum_length stencil linkcentered

electrokinetics 1 density 0.0 D 0.0 valency $valency 0 0

# Set up the charged boundaries 

electrokinetics boundary charge_density [expr $sigma/$agrid] wall normal 0 0 1 d $padding 0 0 direction outside
electrokinetics boundary charge_density [expr $sigma/$agrid] wall normal 0 0 -1 d -[expr $padding+$width] 0 0 direction outside

# Integrate the system

for { set i 0 } { i < 100 } { incr i } {
    part 0 pos [expr $padding + 0.01*$i * $width] 0 0 q 1.0
    integrate 0
    puts [part 0 pr pos f]
}
