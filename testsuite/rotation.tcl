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

proc require_feature {feature} {
    global errf
    if { ! [regexp $feature [code_info]]} {
	set f [open $errf "w"]
	puts $f "not compiled in: $feature"
	close $f
	exit -42
    }
}

puts "----------------------------------------------"
puts "- Testcase rotation.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "----------------------------------------------"
require_feature "ROTATION"
set epsilon 5e-4
setmd temp 0
setmd gamma 0
setmd time_step 0.001
setmd skin 0.5

proc read_data {file} {
    set f [open $file "r"]
    while {![eof $f]} { blockfile $f read auto}
    close $f
}

if { [catch {
    read_data "gb_system.data"

    # to ensure force recalculation
    invalidate_system
    
    inter 0 0 gay-berne 1.0 1.0 4.0 3.0 5.0 2.0 1.0
    
    set GBeng_0 [expr [analyze energy gb 0 0]]
    set toteng_0 [analyze energy total]
    if { [expr abs($toteng_0 - $GBeng_0)] > $epsilon } {
	error "system has unwanted energy contribution, i.e. U_GB != U_total"
    }
  
    integrate 50

    # check the conservation of the total energy
    set toteng [analyze energy total]
    set rel_eng_error [expr abs(($toteng_0 - $toteng)/$toteng)]
    puts "total energy deviation: $rel_eng_error"
    if { $rel_eng_error > $epsilon } {
	error "relative energy error is too large"
    }
    
    # check new GB energy against expected value
    set GB_expected -2971.72
    set GBeng [expr [analyze energy gb 0 0]]
    set rel_eng_error [expr abs(($GBeng - $GB_expected)/$toteng)]
    puts "   GB energy deviation: $rel_eng_error"
    if { $rel_eng_error > $epsilon } {
	error "relative energy error is too large"
    }

} res ] } {
    error_exit $res
}

exec rm -f $errf
exit 0
