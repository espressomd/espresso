# Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#  Max-Planck-Institute for Polymer Research, Theory Group
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

require_feature "SD" 
require_max_nodes_per_side 1
require_cudadevice
require_feature "EXTERNAL_FORCES"
require_feature "SD_NOT_PERIODIC"
require_feature "SD_FF_ONLY"


proc mobi_phi {dist} {
    set a [setmd sd_radius]
    set r [expr $dist/$a]
    set eta [setmd sd_visc]
    # BD:
    #return [expr 1./(6*[PI]*$eta*$a)]
    # SD:
    return [expr (1 - ( 3./4./$r + .5/$r/$r/$r))/(6*[PI]*$eta*$a)]
}

proc mobi_r {dist} {
    set a [setmd sd_radius]
    set r [expr $dist/$a]
    set eta [setmd sd_visc]
    return [expr (1 - ( 3./2./$r - 1./$r/$r/$r))/(6*[PI]*$eta*$a)]
}

set timesteps 1000
setmd time_step [expr [PI]/$timesteps]
setmd sd_visc 1
setmd sd_radius 0.6666666666666666666666
setmd box_l 10 10 10
setmd skin 0.1

part 0 pos  1 0 0
part 1 pos -1 0 0

setmd sd_visc [expr [mobi_phi 2]]

set m2 [mobi_phi 2]
if {$m2 > 1.001 || $m2 < 0.999} {
    puts "failed to set viscosity"
    puts "aborting"
    error_exit
}

thermostat sd 0

for {set i 0} {$i < $timesteps} {incr i} {
    set ppos [part 0 print pos]
    part 0 ext_force        [lindex $ppos 1] [expr -[lindex $ppos 0]] 0
    part 1 ext_force [expr -[lindex $ppos 1]]       [lindex $ppos 0] 0
    integrate_sd 1
}

set pos "[part 0 print pos] [part 1 print pos]" 
# small deviation from $\pm1$ in the positons, because euler is not precise
# radius is about 0.49\% to big with 1000 timesteps and 5.05\% for 100 timesteps 
set target {-1.0049 0 0 1.0049 0 0}
set prec 0.005

for {set d 0} {$d < 6} {incr d} {
    set err [expr abs([lindex $pos $d] - [lindex $target $d])]
    if {$err > $prec } {
	puts "Error to large ($err). Something is wrong with Stokesian Dynamics."
	puts "Positions are: \n$pos\nbut should be:\n$target"
	error_exit
    }
}


part 0 pos  1 0 0 ext_force 1 0 0
part 1 pos -1 0 0 ext_force -1 0 0
set timesteps 3
setmd sd_visc 1
setmd time_step [expr 10./$timesteps]
setmd sd_radius 1
set r 2
for {set i 0} {$i < $timesteps} {incr i} {
    set r [expr $r + 2*[mobi_r $r]*[setmd time_step]]
}
integrate_sd $timesteps
set dist [expr abs( [lindex  [part 0 print pos] 0] - [lindex [part 1 print pos] 0]) ]
set err [expr abs($dist - $r )]
if {$err > 0.00000000001} {
    puts "Test of mobility in radial direction failed."
    puts "Error was $err"
    error_exit
}
