#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
    PLATFORM=`uname -s`; if [ "$1" != "" ]; then NP=$1; else NP=1; fi
# OSF1 \
    if test $PLATFORM = OSF1; then  exec dmpirun -np $NP $TCLMD_SOURCE/$PLATFORM/tcl_md $0 $*
# AIX \
    elif test $PLATFORM = AIX; then exec poe $TCLMD_SOURCE/$PLATFORM/tcl_md $0 $* -procs $NP
# Linux \
    else export EF_ALLOW_MALLOC_0=1; lamboot; exec mpirun -np $NP -nsigs $TCLMD_SOURCE/$PLATFORM/tcl_md $0 $*;
# \
    fi;

#############################################################
#                                                           #
#  Analyze for:                                                         #
#  Test System 3: Polyelectrolyte Solution                  #
#                 (Poor Solvent)                            #
#                                                           #
#  Created:       26.03.2003 by HL                          #
#  Last modified: 26.03.2003 by HL                          #
#                                                           #
#############################################################

puts " "
puts "======================================================="
puts "=                 pe_analyze.tcl                      ="
puts "======================================================="
puts " "

# Simulation to analyze 
# Assumes Configuration files in blockfileformat: $name$ident.####
set name  "pe_solution"
set ident "_t4"

# Specify which configurations to analyze
set anaident ""
set start 60
set end   99
set step  1

set tcl_precision 8

# read in configuration and add up the analyzation values
set config_cnt 0
for { set config $start } { $config <= $end } { incr config $step} {

    # Read in configuration
    set file [open "$name$ident.[format %04d $config]" r]
    while { [blockfile $file read auto] != "eof" } {}
    close $file

    # Do analyzations
    set dist [analyze distribution {2} {0 1} 1.0 300.0 600 1 1]

    # Add up analyzation values ibn result
    if { $config == $start } {
	set result "[lindex $dist 1]"
    } else {
	set temp "[lindex $dist 1]"
	for { set i 0 } {$i < [llength $temp]} {incr i } {
	    set ind [list $i 1]
	    lset result $ind [expr [lindex $result $ind] + [lindex $temp $ind] ]
	}
    }
    puts -nonewline "Analyze ($start-$end) now: $name$ident.[format %04d $config]\r"
    flush stdout
    incr config_cnt
}
# average the analyzation values
for { set i 0 } {$i < [llength $temp]} {incr i } {
    set ind [list $i 1]
    lset result $ind [expr [lindex $result $ind]/$config_cnt]
}

# write results to file
set file [open "$name$ident.iondis$anaident.dat" w]
puts $file "\# r\tP(r)"
for { set i 0 } {$i < [llength $temp]} {incr i } {
    puts $file "[lindex $result $i]"
}
close $file

puts "\nFinished."
