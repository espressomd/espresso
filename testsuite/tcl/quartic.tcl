# Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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

source "tests_common.tcl"

if {[has_feature "LEES_EDWARDS"]} {
    require_max_nodes_per_side 2
}

puts "------------------------------------------"
puts "- Testcase quartic.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "------------------------------------------"

set epsilon 1e-4
set k0 30.
set k1 10.
set r0 2.5
set tolerance 1e-5

set box_l 10.

set n_pairs 300

setmd time_step 0.006
setmd skin      1.0
thermostat off
setmd box_l $box_l $box_l $box_l

proc vecdot {v1 v2} {
    if {[llength $v1]!=[llength $v2]} {
    error "Vectors must be of equal length"
    }
   set dot [expr [lindex $v1 0]*[lindex $v2 0]+[lindex $v1 1]*[lindex $v2 1]+[lindex $v1 2]*[lindex $v2 2]]
   return $dot
}

proc vecadd {c1 v1 c2 v2} {
    if {[llength $v1]!=[llength $v2]} {
	error "Vectors must be of equal length"
    }
    for { set i 0 } { $i < [llength $v1] } { incr i } {
	lappend ret [expr $c1*[lindex $v1 $i] + $c2*[lindex $v2 $i]]
    }
    return $ret
}

inter 0 quartic $k0 $k1 $r0

set energy 0.0
set lastid 0

for { set i 0 } { $i < $n_pairs } { incr i } {
    set posx [expr $box_l*[t_random]]
    set posy [expr $box_l*[t_random]]
    set posz [expr $box_l*[t_random]]
    set r1 [list $posx $posy $posz]

    part $lastid pos $posx $posy $posz
    incr lastid

    set posx [expr $box_l*[t_random]]
    set posy [expr $box_l*[t_random]]
    set posz [expr $box_l*[t_random]]
    set r2 [list $posx $posy $posz]

    part $lastid pos $posx $posy $posz bond 0 [expr $lastid-1]
    incr lastid    

    set r [vecadd 1.0 $r2 -1.0 $r1]
    
    for { set j 0 } { $j < 3 } { incr j } {
	set dx [lindex $r $j] 
	if { $dx >= [expr 0.5*$box_l] } {
	    lset r $j [expr $dx - $box_l]
	}
	if { $dx <= [expr -0.5*$box_l] } {
	    lset r $j [expr $dx + $box_l]
	}
    }
    puts "r $r"

    set dist [expr { sqrt([vecdot $r $r]) }]
    set dr [expr $dist - $r0]

    set energy [expr $energy + 0.5 * $k0 * $dr * $dr + 0.25 * $k1 * $dr*$dr*$dr*$dr]

    lappend drs $dr
    lappend forces [expr abs($k0 * $dr + $k1 * $dr*$dr*$dr)]    
}

integrate 0

set rms 0
if { [catch {
    for { set i 0 } { $i < $n_pairs } { incr i } {
	set r1 [part [expr 2*$i+0] pr pos]
	set r2 [part [expr 2*$i+1] pr pos]
	set r [vecadd 1.0 $r1 -1.0 $r2]
	set r_l [expr { sqrt( [vecdot $r $r] ) }]

	set f1 [part [expr 2*$i+0] pr f]
	set f1_l [expr { sqrt( [vecdot $f1 $f1] ) }]
	set f2 [part [expr 2*$i+1] pr f]
	set f [vecadd 1.0 $f1 1.0 $f2]

	if { [expr { sqrt([vecdot $f $f]) }] > $tolerance } {
	    error "Forces are not oppsite."
	}

	if {  [expr abs([vecdot $r $f1]/($r_l*$f1_l)) - 1] > $tolerance } {
	    error "Force is not parallel to connection vector."
	}

	set ref_force [lindex $forces $i]
	if { [expr abs(($f1_l - $ref_force)/$ref_force)] > $tolerance } {
	    puts "i $i ref_force $ref_force f1_l $f1_l dr [lindex $drs $i]"	    
	    error "Wrong force magnitude."
	}
    }


    set cureng [analyze energy bonded 0 0]
    set energy_error [expr abs(($energy-$cureng)/$energy)]

    if { $energy_error > $tolerance } {
	error "Energy error $energy_error larger than tolerance."
    }
} res ] } {
    error_exit $res
}

exit 0

