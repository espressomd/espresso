# Copyright (C) 2010,2011,2012,2013 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#    Max-Planck-Institute for Polymer Research, Theory Group
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
source "../scripts/bundle.tcl"

require_feature "PARTIAL_PERIODIC"
require_feature "CONSTRAINTS"
require_feature "LENNARD_JONES"

puts "-------------------------------------------------"
puts "- Testcase constraints.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "-------------------------------------------------"

set data_file "constraints_system.data"
set epsilon 1e-4
thermostat off
setmd time_step 0.01
setmd skin 0.05
set energy 0

proc setup_system {file new} {
    set box_l  10.0
    set density 0.3
    set center [expr $box_l/2.0]
    set save   2.0
    global energy

   # simulation box
    setmd box_l $box_l $box_l $box_l
    setmd periodic 0 0 0

    # sphere (type 1)
    set rad           [expr $center-$save]
    set sphere_volume [expr 4.18879*$rad*$rad*$rad]
    constraint sphere center $center $center $center radius $rad type 1
    # cubic box build with walls (type 3)
    constraint wall normal  1  0  0 dist [expr $save] type 3
    constraint wall normal -1  0  0 dist [expr -($box_l-$save)]  type 3
    constraint wall normal  0  1  0 dist [expr $save] type 3
    constraint wall normal  0 -1  0 dist [expr -($box_l-$save)]  type 3
    constraint wall normal  0  0  1 dist [expr $save] type 3
    constraint wall normal  0  0 -1 dist [expr -($box_l-$save)]  type 3
    # cylinder with caps (type 5)
    constraint cylinder center $center $center $center axis 0 0 1 \
	radius $rad length $box_l direction -1 type 5

    # set up three types of particles
    if { $new ==1 } {
	set n_part [expr int($sphere_volume*$density)]
	bundle_counterion_setup $n_part [expr $rad-0.3] 0.0 \
	    "$center $center $center" 0 0
	bundle_counterion_setup $n_part [expr $rad-0.3] 0.0 \
	    "$center $center $center" 2 $n_part
	bundle_counterion_setup $n_part [expr $rad-0.3] 0.0 \
	    "$center $center $center" 4 [expr 2*$n_part]
    }

    # interactions between the particles
    inter 0 0 lennard-jones 1.0 1.0 1.12246 0.25 0.0
    inter 2 2 lennard-jones 1.0 1.0 1.12246 0.25 0.0
    inter 4 4 lennard-jones 1.0 1.0 1.12246 0.25 0.0
    # interactions between particles and constraints
    inter 0 1 lennard-jones 1.0 1.0 1.12246 0.25 0.0
    inter 2 3 lennard-jones 1.0 1.0 1.12246 0.25 0.0
    inter 4 5 lennard-jones 1.0 1.0 1.12246 0.25 0.0
    inter 4 3 lennard-jones 1.0 1.0 1.12246 0.25 0.0

    if { $new ==1 } {
	# equilibrate system with lj-cap potential
	set cap 100
	inter forcecap $cap
	setmd time_step 0.0001
	for { set i 0 } { $i < 100 } { incr i } {
	    integrate 100
	    puts -nonewline "warmup $i from 1000; energy [lindex [lindex [analyze energy] 0] 1]\r"
	    flush stdout
	    set cap [expr $cap+10]
	    inter forcecap $cap
	}
	# turn of capping and write data
	inter forcecap 0
	integrate 0
	puts "\nintegration done."
	set energy [analyze energy]
	puts "Energy without cap: $energy"
	write_data $file
	puts "data written to $file\n"
    }
}

proc read_data {file} {
    global energy
    set f [open $file "r"]
    while {![eof $f]} { blockfile $f read auto}
    close $f
}

proc write_data {file} {
    global energy
    set f [open $file "w"]
    blockfile $f write variable box_l
    blockfile $f write tclvariable energy
    blockfile $f write particles {id type pos q f }
    close $f
}

if { [catch {
    global energy

    if { ![file exists $data_file ] } {
	setup_system "constraints_system.data" 1
	exit
    } else {
	setup_system "constraints_system.data" 0
	read_data $data_file
    }

    for { set i 0 } { $i <= [setmd max_part] } { incr i } {
	set F($i) [part $i pr f]
    }
 
    integrate 0

    set new_energy [analyze energy]
    # check energies
    set maxde 0 
    for { set i 2 } { $i < [llength $new_energy] } { incr i } {
	set de [expr abs([lindex [lindex $energy $i] 3]-[lindex [lindex $new_energy $i] 3]) ]
	if { $de > $maxde } { set maxde $de }
	if { $de > $epsilon } {
	    puts "Energy Error for [lindex $new_energy $i] vs. [lindex $energy $i]: deviation $maxde"
	    #	    error "energy error too large"
	}
    }
    puts "maximal energy deviation $maxde"
    # check forces
    set maxdx 0
    set maxpx 0
    set maxdy 0
    set maxpy 0
    set maxdz 0
    set maxpz 0
    for { set i 0 } { $i <= [setmd max_part] } { incr i } {
	set resF [part $i pr f]
	set tgtF $F($i)
	set dx [expr abs([lindex $resF 0] - [lindex $tgtF 0])]
	set dy [expr abs([lindex $resF 1] - [lindex $tgtF 1])]
	set dz [expr abs([lindex $resF 2] - [lindex $tgtF 2])]

	if { $dx > $maxdx} {
	    set maxdx $dx
	    set maxpx $i
	}
	if { $dy > $maxdy} {
	    set maxdy $dy
	    set maxpy $i
	}
	if { $dz > $maxdz} {
	    set maxdz $dz
	    set maxpz $i
	}
    }
    puts "maximal force deviation in x $maxdx for particle $maxpx, in y $maxdy for particle $maxpy, in z $maxdz for particle $maxpz"
    if { $maxdx > $epsilon || $maxdy > $epsilon || $maxdz > $epsilon } {
	if { $maxdx > $epsilon} {puts "force of particle $maxpx: [part $maxpx pr f] != $F($maxpx)"}
	if { $maxdy > $epsilon} {puts "force of particle $maxpy: [part $maxpy pr f] != $F($maxpy)"}
	if { $maxdz > $epsilon} {puts "force of particle $maxpz: [part $maxpz pr f] != $F($maxpz)"}
	error "force error too large"
    }
} res ] } {
    error_exit $res
}

exit 0
