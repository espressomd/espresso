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
#  Copyright (c) 2002-2005; all rights reserved unless otherwise stated.
# 
puts "-------------------------------------------"
puts "- Testcase el2d.tcl running on [format %02d [setmd n_nodes]] nodes: -"
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
require_feature "PARTIAL_PERIODIC"

set epsilon 1e-3
setmd temp 0
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
    blockfile $f write tclvariable energy
    blockfile $f write interactions
    blockfile $f write particles {id pos q f}
    close $f
}

if { [catch {
    read_data "el2d_system.data"

    for { set i 0 } { $i <= [setmd max_part] } { incr i } {
	set F($i) [part $i pr f]
    }

    # buffer the balanced node_grid for domain_decomposition/P3M
    set balance_ng [setmd node_grid]

    #########################################
    ####      check MMM2D first
    #########################################

    setmd periodic 1 1 0
    cellsystem layered [expr 24/[setmd n_nodes]]
    inter coulomb 1.0 mmm2d 1e-4
    puts -nonewline "MMM2D   node_grid=[setmd node_grid]  "
    flush stdout

    # to ensure force recalculation
    invalidate_system
    integrate 0

    # here you can create the necessary snapshot
    if { 0 } {
	# makes sure that particles are only in the lower 90\% of the box for elc
	set maxz 0
	set minz 1e100
	for {set i 0} {$i <= [setmd max_part]} {incr i} {
	    set pos [part $i pr folded]
	    set pz [lindex $pos 2]
	    if {$pz > $maxz} {set maxz $pz}
	    if {$pz < $minz} {set minz $pz}
	}
	puts "$minz<=z<=$maxz"
	if { $maxz > [expr 0.9*[lindex [setmd box_l] 2]] } {
	    puts "rescaling necessary"
	    # rescale coordinates
	    for {set i 0} {$i <= [setmd max_part]} {incr i} {
		set pos [part $i pr folded]
		set px [lindex $pos 0]
		set py [lindex $pos 1]
		set pz [expr 0.9*[lindex $pos 2]]
		part $i pos $px $py $pz
	    }
	}
	# the high precision requires a lower n_layers
	cellsystem layered 10
	inter coulomb 1.0 mmm2d 1e-20
	integrate 0
	inter coulomb 0.0
	write_data "el2d_system.data"
    }

    ############## RMS force error for MMM2D

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
	error "MMM2D force error too large"
    }

    set cureng [lindex [analyze energy coulomb] 0]
    set toteng [analyze energy total]

    if { [expr abs($toteng - $cureng)] > $epsilon } {
	error "system has unwanted energy contributions of [format %e [expr $toteng - $cureng]]"
    }

    set rel_eng_error [expr abs(($toteng - $energy)/$energy)]
    puts "relative energy deviations: $rel_eng_error"
    if { $rel_eng_error > $epsilon } {
	error "relative energy error too large"
    }

    inter coulomb 0.0

    #########################################
    ####      check P3M + ELC first
    #########################################

    cellsystem domain_decomposition

    eval setmd node_grid $balance_ng
    puts -nonewline "P3M+ELC node_grid=[setmd node_grid]  "
    flush stdout

    setmd periodic 1 1 1
    inter coulomb 1.0 p3m 16.3929265333 32 4 0.142069354322 9.78886014586e-05
    inter coulomb epsilon metallic n_interpol 32768 mesh_off 0.5 0.5 0.5
    inter coulomb elc 1e-4 [expr 0.1*[lindex [setmd box_l] 2]]

    # to ensure force recalculation
    invalidate_system
    integrate 0

    ############## RMS force error for P3M + ELC

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
	error "P3M+ELC force error too large"
    }

    set cureng [lindex [analyze energy coulomb] 0]
    set toteng [analyze energy total]

    puts [analyze energy]

    if { [expr abs($toteng - $cureng)] > $epsilon } {
	error "system has unwanted energy contributions of [format %e [expr $toteng - $cureng]]"
    }

    set rel_eng_error [expr abs(($toteng - $energy)/$energy)]
    puts "relative energy deviations: $rel_eng_error"
    if { $rel_eng_error > $epsilon } {
	error "relative energy error too large"
    }
} res ] } {
    error_exit $res
}

exec rm -f $errf
exit 0
