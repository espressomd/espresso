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
#  Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
# 
set errf [lindex $argv 1]

proc error_exit {error} {
    global errf
    set f [open $errf "w"]
    puts $f "Error occured: $error"
    close $f
    exit -666
}

puts "------------------------------------------------"
puts "- Testcase thermostat.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "------------------------------------------------"
set epsilon 3e-2
setmd temp 1
setmd gamma 1
setmd time_step 0.01
setmd skin 0.5
set n_part 100
set maxstep 100
set intstep [expr round(50*$maxstep/100)]

proc read_data {file} {
    set f [open $file "r"]
    while {![eof $f]} { blockfile $f read auto}
    close $f
}

proc write_data {file} {
    set f [open $file "w"]
    blockfile $f write variable box_l
    blockfile $f write particles {id pos v omega} 
    close $f
}
if { [catch {
    if { [regexp "ROTATION" [code_info]] } {
	puts "rotation found, 6 degrees of freedom"
	set deg_free 6
	set filename "thermostat_rot.data"
    } else {
	puts "no rotation found, 3 degrees of freedom"
	set deg_free 3
	set filename "thermostat.data"
    }
    read_data $filename

    set eng0    [analyze energy kin]
    set temp0   [expr $eng0/$n_part/($deg_free/2.)]
    set curtemp1 0

    for {set i 0} { $i < 100} { incr i } {
    integrate $intstep
	puts -nonewline " $i percent done\r"
	flush stdout
	set toteng [analyze energy total]
	set cureng [analyze energy kin] 
	set curtemp [expr $cureng/$n_part/($deg_free/2.)] 

	if { [expr abs($toteng - $cureng)] > $epsilon } {
	    error "system has unwanted energy contributions"
	}
	set curtemp1 [expr $curtemp1 + $curtemp]
    }
    set curtemp1 [expr $curtemp1/$maxstep]
    # here you can create a new snapshot
    # write_data $filename


    set rel_temp_error [expr abs(( [setmd temp] - $curtemp1)/$curtemp1)]
    puts "thermostat temperature:          [setmd temp]"
    puts "measured temperature:            $curtemp1"
    puts "relative temperature deviations: $rel_temp_error"
    if { $rel_temp_error > $epsilon } {
	error "relative temperature error too large"
    }
  
} res ] } {
    error_exit $res
}

exec rm -f $errf
exit 0