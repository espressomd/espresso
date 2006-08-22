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
puts "-------------------------------------------"
puts "- Testcase p3m.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "-------------------------------------------"
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

require_feature "ELECTROSTATICS"

set epsilon 1e-3
thermostat off
setmd time_step 0.01
setmd skin 0.05


proc read_data {file} {
    set f [open $file "r"]
    while {![eof $f]} { blockfile $f read auto}
    close $f
}

proc write_data {file} {
    set f [open $file "w"]
    blockfile $f write variable box_l
    blockfile $f write tclvariable {energy pressure}
    blockfile $f write interactions
    blockfile $f write particles {id pos q f}
    close $f
}

if { [catch {
    read_data "p3m_system.data"

    for { set i 0 } { $i <= [setmd max_part] } { incr i } {
	set F($i) [part $i pr f]
    }
    ############## P3M-specific part
    # the P3M parameters are stored in p3m_system.data

    # to ensure force recalculation
    invalidate_system
    integrate 0

    # here you can create the necessary snapshot
    if { 0 } {
	inter coulomb 1.0 p3m tune accuracy 1e-4
	integrate 0

	write_data "p3m_system.data"
    }

    ############## end

    puts [analyze energy]
    puts [analyze pressure]

    set cureng [lindex [analyze   energy coulomb] 0]
    set curprs [lindex [analyze pressure coulomb] 0]

    set rel_eng_error [expr abs(($cureng - $energy)/$energy)]
    puts "relative energy deviations: $rel_eng_error"
    if { $rel_eng_error > $epsilon } {
	error "relative energy error too large"
    }

    set rel_prs_error [expr abs(($curprs - $pressure)/$pressure)]
    puts "relative pressure deviations: $rel_prs_error"
    if { $rel_prs_error > $epsilon } {
	error "relative pressure error too large"
    }

    ############## end, here RMS force error for P3M

    set rmsf 0
    for { set i 0 } { $i <= [setmd max_part] } { incr i } {
	set resF [part $i pr f]
	set tgtF $F($i)
	set dx [expr abs([lindex $resF 0] - [lindex $tgtF 0])]
	set dy [expr abs([lindex $resF 1] - [lindex $tgtF 1])]
	set dz [expr abs([lindex $resF 2] - [lindex $tgtF 2])]

	set rmsf [expr $rmsf + $dx*$dx + $dy*$dy + $dz*$dz]
    }
    set rmsf [expr sqrt($rmsf/[setmd n_part])]
    puts "rms force deviation $rmsf"
    if { $rmsf > $epsilon } {
	error "force error too large"
    }
} res ] } {
    error_exit $res
}

exec rm -f $errf
exit 0
