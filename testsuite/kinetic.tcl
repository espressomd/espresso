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

puts "---------------------------------------------"
puts "- Testcase kinetic.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "---------------------------------------------"
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
require_feature "PARTIAL_PERIODIC"
require_feature "EXTERNAL_FORCES"

set epsilon 1e-8
setmd temp 0
setmd gamma 0
setmd time_step 1
setmd skin 0
set   box_l 200
setmd box_l $box_l $box_l $box_l
setmd periodic 0 0 0

set bjerrum 10.0
set maxr [expr floor(0.5*$box_l)-1]
set q1 +1.0
set q2 +1.0
set x2 +5.0
set v2 -1.0

set maxstep 50
set maxx 0
set maxf 0
set maxe 0
set maxk 0


if { [catch {
    inter coulomb $bjerrum dh 0 $maxr
    part 0 pos 0.0 0.0 0.0 type 1 q $q1 fix
    part 1 pos $x2 0.0 0.0 type 1 q $q2 v $v2 0.0 0.0

    for {set i 0} {$i<$maxstep} {incr i} {
	integrate 10

	set pos [part 1 print pos]; set x [lindex $pos 0]; set y [lindex $pos 1]; set z [lindex $pos 2]
	if { [expr abs($y)] > $epsilon || [expr abs($z)] > $epsilon } {
	    error "unexpected shear occured ($y / $z)"
	} elseif { $x > $maxr } { break }
	set vel [part 1 print v]; set vx [lindex $vel 0]; set vy [lindex $vel 1]; set vz [lindex $vel 2]
	if { [expr abs($vy)] > $epsilon || [expr abs($vz)] > $epsilon } {
	    error "unexpected shear velocities occured ($vy / $vz)"
	}
	set curf [part 1 print f]; set curfx [lindex $curf 0]; set curfy [lindex $curf 1]; set curfz [lindex $curf 2]
	if { [expr abs($curfy)] > $epsilon || [expr abs($curfz)] > $epsilon } {
	    error "unexpected shear forces occured ($curfy / $curfz)"
	}
	set curedh [analyze energy coulom]; set curekin [analyze energy kinetic]; set curet [analyze energy total]
	if { [expr abs($curet - $curedh - $curekin)] > $epsilon } {
	    error "system has unwanted energy contributions of [format %e [expr $curet - $curedh - $curekin]]"
	}

	set thx [analyze mindist]; set rel_error [expr abs(($x - $thx)/$thx)]
	if { $rel_error > $epsilon } {
	    error "relative error $rel_error in the minimal distance too large ($thx / $x)"
	} elseif { $rel_error > $maxx } { set maxx $rel_error }
	set thfx [expr $bjerrum*$q1*$q2/($x*$x)]; set rel_error [expr abs(($curfx - $thfx)/$thfx)]
	if { $x > $maxr } { set thfx 0.0; set rel_error [expr abs($curfx - $thfx)] }
	if { $rel_error > $epsilon } {
	    error "relative error $rel_error in the coulomb force too large ($thfx / $curfx)"
	} elseif { $rel_error > $maxf } { set maxf $rel_error }
	set thedh [expr $bjerrum*$q1*$q2/$x]; set rel_error [expr abs(($curedh - $thedh)/$thedh)]
	if { $x > $maxr } { set thedh 0.0; set rel_error [expr abs($curedh - $thedh)] }
	if { $rel_error > $epsilon } {
	    error "relative error $rel_error in the coulomb energy too large ($thedh / $curedh)"
	} elseif { $rel_error > $maxe } { set maxe $rel_error }
	set thekin [expr 0.5*$vx*$vx]; set rel_error [expr abs(($curekin - $thekin)/$thekin)]
	if { $rel_error > $epsilon } {
	    error "relative error $rel_error in the kinetic energy too large ($thek / $curekin)"
	} elseif { $rel_error > $maxk } { set maxk $rel_error }
    }
    puts "maximal deviation in the minimal distance = $maxx, in the force = $maxf, in the coulomb energy = $maxe, in the kinetic energy = $maxk"
} res ] } {
    error_exit $res
}

exec rm -f $errf
exit 0
