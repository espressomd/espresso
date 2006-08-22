#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
    exec $ESPRESSO_SOURCE/Espresso $0 $*
# 
#  This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
#  It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
#  and by which you are legally bound while utilizing this file in any form or way.
#  There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  You should have received a copy of that license along with this program;
#  if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
#  write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
#  Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
# 
set errf [lindex $argv 1]

proc error_exit {error} {
    global errf
    set f [open $errf "w"]
    puts $f "Error occured: $error"
    close $f
    exit -666
}

proc require_feature {feature} {
    global errf
    if { ! [regexp $feature [code_info]]} {
	set f [open $errf "w"]
	puts $f "not compiled in: $feature"
	close $f
	exit -42
    }
}

puts "----------------------------------------"
puts "- Testcase comfixed.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "----------------------------------------"

require_feature "COMFIXED"
require_feature "PARTIAL_PERIODIC"

set epsilon 1e-4
thermostat off
setmd time_step 1
setmd skin 0

if { [setmd n_nodes] != 1 } {
    puts "Testcase comfixed.tcl does not run on more than one node"
    exec rm -f $errf
    exit 0
}


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
    # to ensure force recalculation
    invalidate_system

    ############## comfixed-specific part

    setmd box_l 10. 10. 10.
    setmd skin  0.5
    setmd time_step 0.005
    setmd periodic  0 0 0
    thermostat langevin 1.0 1.0

    inter 0 0 lennard-jones 1.0 1.0 1.12246 0.25 0.0
    constraint sphere center 5. 5. 5.  radius 4. type 5
    inter 0 5 lennard-jones 1.0 1.0 1.12246 0.25 0 0
    inter 0 0 comfixed 1

    part 0 pos 5 5 5 type 0
    part 1 pos 6 5 5 type 0
    part 2 pos 6 6 4 type 0
    part 3 pos 4 5 7 type 0

    set com_ini [analyze centermass 0]
    integrate 100
    set com_fin [analyze centermass 0]
    set com_dist [bond_length $com_ini $com_fin]

    if { $com_dist > $epsilon } {
	error "center of mass of particles has moved $com_dist"
    }

} res ] } {
    error_exit $res
}

exec rm -f $errf
exit 0
