#!/bin/sh
# tricking... the line after a these comments are interpreted as standard shell script \
    PLATFORM=`uname -s`; if [ "$1" != "" ]; then NP=$1; else NP=2; fi
# OSF1 \
    if test $PLATFORM = OSF1; then  exec dmpirun -np $NP $ESPRESSO_SOURCE/$PLATFORM/Espresso $0 $*
# AIX \
    elif test $PLATFORM = AIX; then exec poe $ESPRESSO_SOURCE/$PLATFORM/Espresso $0 $* -procs $NP
# Linux \
    else export EF_ALLOW_MALLOC_0=1; exec mpirun -np $NP -nsigs $ESPRESSO_SOURCE/$PLATFORM/Espresso $0 $*;
# \
    fi;
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

puts "----------------------------------------"
puts "- Testcase npt.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "----------------------------------------"


set epsilon 1e-4
setmd temp 0
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
    blockfile $f write particles {id type pos f q}
    close $f
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

if { [catch {

require_feature "NPT"
require_feature "LENNARD_JONES"

read_data "npt_lj_system.data"

    ############## lj-specific part



# ----------- Integration Parameters before warmup -----------#
setmd time_step 0.01
setmd skin      0.4
thermostat set langevin 1.0 1.0 
setmd temp      1.0

set p_ext 2.0
set piston_mass 0.0001
set gamma_0 0.5
set gamma_v 0.001
set blen [lindex [setmd box_l] 0]

inter 0 0 lennard-jones 1.0 1.0 1.12246 0.25 0.0

integrate set npt_isotropic $p_ext $piston_mass 
#-cubic_box
thermostat set npt_isotropic 1.0 $gamma_0 $gamma_v 

set avp 0
set avpi 0
set avbl 0
set avbli 0

puts -nonewline "Running NPT Test for uncharged lj system "
for { set t 0 } { $t < 1000 } { incr t } {
    integrate 10
    set totprs [analyze pressure total]
    set avp [expr $avp + $totprs]
    incr avpi

    set avbl [expr $avbl + [lindex [setmd box_l] 0] ]
    incr avbli

    if { $t%100==0 } {
	puts -nonewline ". "
	flush stdout
    }

}
if { [expr abs($avp/(1.0*$avpi) - 2.0) ] > 0.2 } {
    error "ERROR: Average pressure <P> = [expr $avp/(1.0*$avpi)] deviates from imposed pressure P = 2.0"
} else { puts "done" }  

# here you can create the necessary snapshot
#    write_data "npt_lj_system.data"

# ---------------------------------------------------------------------------#
# ---------------------------------------------------------------------------#
# Test NPT with ELECTROSTATICS  #

# Somebody who understands this please feel free to write it #

} res ] } {
    error_exit $res
}

exec rm -f $errf
exit 0



