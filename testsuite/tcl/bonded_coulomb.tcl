# Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#   Max-Planck-Institute for Polymer Research, Theory Group
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

source "tests_common.tcl"

require_feature "ELECTROSTATICS"
require_feature "COULOMB_DEBYE_HUECKEL" off

puts "------------------------------------------"
puts "- Testcase bonded_coulomb.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "------------------------------------------"

set tolerance 1e-5
set bjerrum 2.0
set box_l 10
set n_parts 10
set max_charge 5.0

setmd time_step 0.006
setmd skin      1.0
thermostat off
setmd box_l $box_l $box_l $box_l
cellsystem nsquare

set lastid 0

# Set up randome charge system

for { set i 0 } { $i < $n_parts } { incr i } {
    set posx [expr $box_l*[t_random]]
    set posy [expr $box_l*[t_random]]
    set posz [expr $box_l*[t_random]]

    part $lastid pos $posx $posy $posz q [expr (0.5 - [t_random])*$max_charge]

    incr lastid
}

# Set up bonds

inter 0 bonded_coulomb $bjerrum

for { set i 0 } { $i < $n_parts } { incr i } {
    for { set j [expr $i + 1] } { $j < $n_parts } { incr j } {
	part $i bond 0 $j
    }
}

integrate 0

for { set i 0 } { $i < $n_parts } { incr i } {
    lappend forces_bonded_coulomb [part $i pr f]
}

set energy_bonded_coulomb [analyze energy bonded 0] 

#Disable bonded IA
inter 0 bonded_coulomb 0.0

inter coulomb $bjerrum dh 0.0 $box_l
set energy_dh [analyze energy coulomb]

integrate 0 recalc_forces

set rms 0
for { set i 0 } { $i < $n_parts } { incr i } {
    set f_dh [part $i pr f]
    set f_bonded_coulomb [lindex $forces_bonded_coulomb $i]
    puts "dh $f_dh"
    puts "bc $f_bonded_coulomb"
    set rms [expr $rms + pow([lindex $f_dh 0] - ([lindex $f_bonded_coulomb 0]), 2)]
    set rms [expr $rms + pow([lindex $f_dh 1] - ([lindex $f_bonded_coulomb 1]), 2)]
    set rms [expr $rms + pow([lindex $f_dh 2] - ([lindex $f_bonded_coulomb 2]), 2)]
}

set rms [expr { sqrt($rms/$n_parts) }]

if { [catch {
    if { $rms >= $tolerance } {
	error "Force rms $rms is large than tolerance $tolerance."
    }
    if { [expr abs(($energy_bonded_coulomb - $energy_dh)/$energy_bonded_coulomb)] >= $tolerance } {
	error "Energy error is too large."
    }
} res ] } {
    error_exit $res
}

exit 0
