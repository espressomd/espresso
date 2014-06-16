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

require_feature "TABULATED"

puts "----------------------------------------"
puts "- Testcase tabulated.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "----------------------------------------"

set epsilon 1e-6
thermostat off
setmd time_step 1
setmd skin 0

proc read_data {file} {
    set f [open $file "r"]
    while {![eof $f]} { blockfile $f read auto}
    close $f
}

proc write_data {file} {
    global energy pressure
    set f [open $file "w"]
    set energy [analyze energy total]
    set pressure [analyze pressure total]
    blockfile $f write tclvariable {energy pressure}
    blockfile $f write variable box_l
    blockfile $f write particles {id type pos f}
    close $f
}

if { [catch {
    read_data "tabulated_system.data" 


    for { set i 0 } { $i <= [setmd max_part] } { incr i } {
	set F($i) [part $i pr f]
    }

    ############## lj-specific part
    inter 0 0 tabulated "lj1.tab"
    inter 1 1 tabulated "lj2.tab"
    inter 0 1 tabulated "lj3.tab"

    inter forcecap 1000000000

# ------------------------------------------------#
# These are the lennard jones potentials that are 
# tabulated in the above files 
#    inter 0 0 lennard-jones 1.0 1.0 1.12246 0.25 0.0
#    inter 1 1 lennard-jones 1.3 0.5 2 0.0 0.0
#    inter 0 1 lennard-jones 2.2 1.0 1.12246 0.0 0.5

    integrate 0



    # here you can create the necessary snapshot
    # write_data "tabulated_system.data"

    # ensures that no other forces are on
    set cureng [expr [analyze energy nonbonded 0 0] + [analyze energy nonbonded 0 1] + [analyze energy nonbonded 1 1]]
    # tbrs
    set curprs [expr [lindex [analyze pressure nonbonded 0 0] 0] + \
		[lindex [analyze pressure nonbonded 0 1] 0] + \
		[lindex [analyze pressure nonbonded 1 1] 0]]

    ############## end

#    puts $cureng

    set toteng [analyze energy total]
    set totprs [analyze pressure total]




    if { [expr abs($toteng - $cureng)] > $epsilon } {
	error "system has unwanted energy contributions of [format %e [expr $toteng - $cureng]]"
    }
    if { [expr abs($totprs - $curprs)] > $epsilon } {
	error "system has unwanted pressure contributions of [format %e [expr $totprs - $curprs]]"
    }

    set rel_eng_error [expr abs(($toteng - $energy)/$energy)]
    puts "relative energy deviations: $rel_eng_error"
    if { $rel_eng_error > $epsilon } {
	error "relative energy error too large"
    }

    set rel_prs_error [expr abs(($totprs - $pressure)/$pressure)]
    puts "relative pressure deviations: $rel_prs_error"
    if { $rel_prs_error > $epsilon } {
	error "relative pressure error too large"
    }

    set maxdx 0
    set maxpx 0
    set maxdy 0
    set maxpy 0
    set maxdz 0
    set maxpz 0
    set maxForce 0
    for { set i 0 } { $i <= [setmd max_part] } { incr i } {
	set resF [part $i print force]

	set tgtF $F($i)

	set dx [expr abs([lindex $resF 0] - [lindex $tgtF 0])]
	set dy [expr abs([lindex $resF 1] - [lindex $tgtF 1])]
	set dz [expr abs([lindex $resF 2] - [lindex $tgtF 2])]

	set force [veclen $resF]
	if {$force > $maxForce} {
		set maxForce $force
	}

	if { $dx > $maxdx} {
	    set maxdx $dx
	    set maxpx $i
#	    puts "here [part $i print force]"
#	    puts $tgtF
#	    flush stdout
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
    set maxForce [expr sqrt($maxForce)]

    puts "maximal force deviation in x $maxdx for particle $maxpx, in y $maxdy for particle $maxpy, in z $maxdz for particle $maxpz"
    if { $maxdx/$maxForce > $epsilon || $maxdy/$maxForce > $epsilon || $maxdz/$maxForce > $epsilon } {
	if { $maxdx > $epsilon} {puts "force of particle $maxpx: [part $maxpx pr f] != $F($maxpx)"}
	if { $maxdy > $epsilon} {puts "force of particle $maxpy: [part $maxpy pr f] != $F($maxpy)"}
	if { $maxdz > $epsilon} {puts "force of particle $maxpz: [part $maxpz pr f] != $F($maxpz)"}
	error "force error too large"
    }
} res ] } {
    error_exit $res
}

exit 0
