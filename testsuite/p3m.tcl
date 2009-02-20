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
puts "--------------------------------------------------------------------------------------------------------"
puts "- Testcase p3m.tcl for CHARGES AND/OR  MAGNETIC DIPOLES running on [format %02d [setmd n_nodes]] nodes: -"
puts "--------------------------------------------------------------------------------------------------------"
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


#============================================================

#if Electrostatics is on, check that P3M for charge-charge interactions works fine
if { [regexp "ELECTROSTATICS" [code_info]]} { 
    source "p3m-charges.tcl"
} else  {
    puts "NO ELECTROSTATICS"  
}


#if Magnetostatics is on, check that P3M for magnetic dipole-dipole interactions works fine
if { [regexp "MAGNETOSTATICS" [code_info]]} { 
    if {[setmd n_nodes] == 1} {
	source "p3m-magnetostatics.tcl"
    } {
	puts "magnetostatics only runs on 1 CPU currently"
    }
} else {
    puts "NO MAGNETOSTATICS "  
}


#CROSS-CHECK:  if both magnetostatic and electrostatic are on, check that one does not inerfer the other
#              it coould happen that due to the share of functions, etc, some error was introduced that
#              only shows when both are active at the same time.

#if { [regexp "ELECTROSTATICS" [code_info]]   &&  [regexp "MAGNETOSTATICS" [code_info]] } { 
  # source "p3m-cross-check-charges-dipoles.tcl"
#}


exec rm -f $errf
exit 0
