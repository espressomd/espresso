# Copyright (C) 2010,2011,2012,2013 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#   Max-Planck-Institute for Polymer Research, Theory Group
#  
# This file is part of ESPResSo.
#  
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 

# 
source "tests_common.tcl"

require_feature "LENNARD_JONES"

puts "----------------------------------------"
puts "- Testcase intpbc.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "----------------------------------------"

set epsilon 1e-4
thermostat off
setmd time_step 0.001
setmd skin 0.05

proc read_data {file} {
    set f [open "|gzip -cd $file" "r"]
    while {![eof $f]} { blockfile $f read auto}
    if { [catch { close $f } fid] } { puts "Error while closing $file caught: $fid." }
}

proc write_data {file} {
    global energy pressure verlet_reuse op
    set f [open "|gzip > $file" "w"]
    set energy [analyze energy total]
    set pressure [analyze pressure total];
    set verlet_reuse [setmd verlet_reuse]
    blockfile $f write tclvariable {energy pressure verlet_reuse}
    blockfile $f write variable box_l
    # particle block by hand as we need the OLD positions
    puts $f "{particles {id pos f}"
    for {set i 0} {$i <= [setmd max_part]} {incr i} {
	puts $f "\t\{[part $i pr id] $op($i) [part $i pr f]\}"
    }
    puts $f "}"

    blockfile $f write bonds
    close $f
}

if { [catch {
    ############## integ-specific part
    setmd box_l     99 99 99
    inter 0 0 lennard-jones 1.0 1.0 1.12246 0.25 0.0

    set fene_k      30.0
    set fene_r      1.5
    inter 0 fene $fene_k $fene_r

    read_data "intpbc_system.data.gz"

    for { set i 0 } { $i <= [setmd max_part] } { incr i } {
	set F($i) [part $i pr f]
    }

    # copy original positions for writing
    for { set i 0 } { $i <= [setmd max_part] } { incr i } {
    	set op($i) [part $i pr pos]
    }

    integrate 100

    # here you can create the necessary snapshot
    # write_data "intpbc_system.data.gz"

    ############## end

    set toteng [analyze energy total]
    set totprs [analyze pressure total]

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

    puts "verlet reuse is [setmd verlet_reuse], should be $verlet_reuse"
    if { [expr abs([setmd verlet_reuse] - $verlet_reuse)] > $epsilon } {
	error "verlet reuse frequency differs."
    }

} res ] } {
    error_exit $res
}

exit 0
