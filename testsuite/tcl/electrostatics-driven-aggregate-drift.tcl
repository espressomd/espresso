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
require_feature "DIPOLES"
require_feature "PARTIAL_PERIODIC"
require_feature "CONSTRAINTS"
require_feature "ELECTROSTATICS"

puts "------------------------------------------------"
puts "- Testcase electrostatics-driven-aggregate-drift.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "------------------------------------------------"
puts "------------------------------------------------"

set PI 3.141592653589793238
# make real random draw
set cmd "t_random seed"
for {set i 0} {$i < [setmd n_nodes]} { incr i } {
  lappend cmd [expr [pid] + $i] }
eval $cmd

# Two aggregates attraction in an external strong magnetic field
# --------------------------------------------------------------

# Number of particles, should be an even number:
set box 1000
setmd box_l $box $box $box
setmd time_step 0.007
cellsystem layered [expr 25/[setmd n_nodes]]

set n 9
set gamma0 7.984
set kT 2.64
setmd periodic 1 1 0
setmd skin 0
thermostat langevin $kT $gamma0
#cellsystem nsquare
#-no_verlet_lists
set J "4.384 4.384 4.384"
set mass 3.847
set B_mag 10
set q 1
set sigma [expr 5/(2.0*$PI)]

# This structure is supposed to be stabilized towards the hexagonal
# closed packed structure in an external strong magnetic field
# like two linear chains.
# The center of mass is zc = $box/2.0 + 7/3.
part 0 pos [expr $box/2.0] [expr $box/2.0] [expr $box/2.0]
part 1 pos [expr $box/2.0] [expr $box/2.0] [expr $box/2.0+2]
part 2 pos [expr $box/2.0] [expr $box/2.0] [expr $box/2.0+4]
part 3 pos [expr $box/2.0] [expr $box/2.0+2] [expr $box/2.0+1]
part 4 pos [expr $box/2.0] [expr $box/2.0+2] [expr $box/2.0+3]
part 5 pos [expr $box/2.0] [expr $box/2.0+2] [expr $box/2.0+5]
part 6 pos [expr $box/2.0] [expr $box/2.0+4] [expr $box/2.0]
part 7 pos [expr $box/2.0] [expr $box/2.0+4] [expr $box/2.0+2]
part 8 pos [expr $box/2.0] [expr $box/2.0+4] [expr $box/2.0+4]

for {set p 0} {$p<$n} {incr p} {
    set theta [expr $PI*[t_random]]
    set phi [expr 2.0*$PI*[t_random]]
    set dx 0
    set dy 0
    set dz 1
    part $p type 0 dip $dx $dy $dz mass $mass q $q rinertia [lindex $J 0] [lindex $J 1] [lindex $J 2] v 0 0 0 omega_body 0 0 0
}

#Interactions
inter magnetic [expr 20/$kT] dawaanr
inter coulomb [expr 1.0/$kT] mmm2d 1E-6

constraint ext_magn_field 0 0 $B_mag
constraint wall normal 0 0 -1 dist [expr -0.97 * $box] type 1 penetrable 0 reflecting 1
constraint plate height [expr 1.1*$box] sigma -$sigma

set sig [expr 5.0/6.0]; set cut [expr 1.4*1.12246*$sig]
set eps 5; set shift [expr 0.25*$eps]
inter 0 0 lennard-jones $eps $sig $cut $shift 0

set sig [expr 5.0/6.0]; set cut [expr 0.4*1.12246*$sig]
set eps 5; set shift [expr 0.25*$eps]
inter 0 1 lennard-jones $eps $sig $cut $shift 0

#setmd time 0
# Dipoles should align along the external magnetic field
#integrate 200
# Now they could start a translational motion
#for {set p 0} {$p<$n} {incr p} {
#    part $p unfix
#}

set q_tot [expr $q*$n]
set F [expr 2*$PI*$sigma*$q_tot]
set gamma_tot [expr $gamma0*$n]
set mass_tot [expr $mass*$n]

prepare_vmd_connection "vmdfile" 10000

# Let's start from zero again
setmd time 0

for {set i 0} { $i < 6000 } { incr i} {
    integrate 10
    imd positions
}

# The aggregate macroscopic parameters
# Center of mass Z component
set pos_c_z 0
# Total magnetic moment along the Z axis
#set mag_total_z 0
for {set p 0} {$p < $n} {incr p} {
    set pos [part [expr $p] print pos]
    set dip [part [expr $p] print dip]
    set pos_c_z [expr $pos_c_z + [lindex $pos 2]]
    #set mag_total_z [expr $mag_total_z + [lindex $dip 2]]
}

set t [setmd time]
set pos_c_z [expr $pos_c_z / $n - $box / 2.0]
set pos_c_z_exp [expr ($F/pow($gamma_tot,2))*($mass_tot*(exp(-$gamma_tot*$t/$mass_tot)-1)+$gamma_tot*$t)+7/3.0];

if { abs($pos_c_z - $pos_c_z_exp) / $pos_c_z_exp > 0.2 } {
    error_exit "Z-drift deviation is too large: $pos_c_z vs the expected value $pos_c_z_exp"
}

exit 0
