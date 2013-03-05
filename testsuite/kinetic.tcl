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

require_feature "ADRESS" off
require_feature "ELECTROSTATICS"
require_feature "PARTIAL_PERIODIC"
require_feature "EXTERNAL_FORCES"

puts "---------------------------------------------"
puts "- Testcase kinetic.tcl running on [format %02d [setmd n_nodes]] nodes: -"
puts "---------------------------------------------"

set epsilon 1e-8
thermostat off
setmd time_step 1
setmd skin 0

set   bjerrum  10.0
set   maxr  100
set   box_l [expr 2*[setmd n_nodes]*($maxr+1)]
setmd box_l $box_l $box_l $box_l
setmd periodic 0 0 0

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
	integrate 1

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
	set curedh [analyze energy coulomb]; set curekin [analyze energy kinetic]; set curet [analyze energy total]
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
	    error "relative error $rel_error in the kinetic energy too large ($thekin / $curekin)"
	} elseif { $rel_error > $maxk } { set maxk $rel_error }
    }
    puts "maximal deviation in the minimal distance = $maxx, in the force = $maxf, in the coulomb energy = $maxe, in the kinetic energy = $maxk"
} res ] } {
    error_exit $res
}

exit 0
