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
#  Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
# 
puts "------------------------------------------"
puts "- Testcase fene.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "------------------------------------------"
setmd time_step 0.006
setmd skin      0.0
setmd gamma     0.0
setmd temp      0.0

set epsilon 1e-4

set errf [lindex $argv 1]
proc error_exit {error} {
    global errf
    set f [open $errf "w"]
    puts $f "Error occured: $error"
    close $f
    exit -666
}

proc read_data {file} {
    if { [string compare [lindex [split $file "."] end] "gz"]==0 } { set f [open "|gzip -cd $file" r]
    } else { set f [open "$file" "r"] }
    while {![eof $f]} { blockfile $f read auto}
    close $f
}

proc write_data {file} {
    global energy pressure
    if { [string compare [lindex [split $file "."] end] "gz"]==0 } { set f [open "|gzip -c - >$file" w]
    } else { set f [open $file "w"] }
    set energy [analyze energy total]
    set pressure [analyze pressure total]
    blockfile $f write tclvariable {energy pressure}
    blockfile $f write variable box_l
    blockfile $f write particles {id pos f}
    blockfile $f write bonds all
    close $f
}

if { [catch {
    ############## fene-specific part
    set fene_k      30.0
    set fene_r      1.5
    setmd box_l     99 99 99
    inter 0 fene $fene_k $fene_r

    read_data "fene_system.data.gz"

    for { set i 0 } { $i <= [setmd max_part] } { incr i } {
	set F($i) [part $i pr f]
    }
    # to ensure force recalculation
    invalidate_system
    integrate 0

    # here you can create the necessary snapshot
    # write_data "fene_system.data2"

    # ensures that no other forces are on
    set cureng [analyze energy fene 0 0]
    # tbrs
    set curprs [lindex [analyze pressure fene 0 0] 0]
    ############## end

    set toteng [analyze energy total]
    set totprs [analyze pressure total]

    if { [expr abs($toteng - $cureng)] > $epsilon } {
	error "system has unwanted energy contributions"
    }
    if { [expr abs($totprs - $curprs)] > $epsilon } {
	error "system has unwanted pressure contributions"
    }

    set rel_eng_error [expr abs(($toteng - $energy)/$energy)]
    puts "relative energy deviations: $rel_eng_error  ($toteng / $energy)"
    if { $rel_eng_error > $epsilon } {
	error "relative energy error too large"
    }

    set rel_prs_error [expr abs(($totprs - $pressure)/$pressure)]
    puts "relative pressure deviations: $rel_prs_error  ($totprs / $pressure)"
    if { $rel_prs_error > $epsilon } {
	error "relative pressure error too large"
    }

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

exec rm -f $errf
exit 0