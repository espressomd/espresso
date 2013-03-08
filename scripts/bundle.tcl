#############################################################
#                                                           #
# bundle.tcl                                                #
# ==========                                                #
#                                                           #
# Functions for bundle simulations.                         #
#                                                           #
#############################################################
#
# Copyright (C) 2010,2012,2013 The ESPResSo project
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
# setup for backbone of cylindrical micelle in sphere 
proc bundle_backbone_setup { n_poly length cyl_rad c_pat h_dist {center "0 0 0"} {types "0 2"} {start "0"} } {
    set mypi          3.141592653589793
    set off  [expr -($h_dist/2.0-0.5)]
    set alpha [expr -(2.0*$mypi/$n_poly)]
    
    for {set i 0} { $i < $n_poly } {incr i} {
	set posx [expr [lindex $center 0] - $length/2.0 + $off + ($i%$h_dist)]
	set posy [expr [lindex $center 1] + $cyl_rad*sin($i*$alpha)]
	set posz [expr [lindex $center 2] + $cyl_rad*cos($i*$alpha)]
	for {set j 0} { $j < $length } {incr j} {
	    set pat_num [expr $j%[llength $c_pat]]
	    if { [lindex $c_pat $pat_num] != 0 } {
		set type [lindex $types 1]
		set charge [lindex $c_pat $pat_num]
	    } else {
		set type [lindex $types 0]
		set charge 0.0
	    }
	    part $start pos $posx $posy $posz type $type q $charge
	    # bonds
	    if { $j > 0 } {
		part [expr $start-1] bond 0 $start
	    }
	    if { $j > 1 } {
		part [expr $start-1] bond 1 [expr $start-2] [expr $start]
	    }


	    set posx [expr $posx+1.0]
	    incr start
	}
    }
}

proc bundle_counterion_setup { n_ci sphere_rad valency {center "0 0 0"} {type "1"} {start "0"} } {
    set s_rad2 [expr $sphere_rad*$sphere_rad]

    for {set i $start} { $i < [expr $n_ci+$start] } {incr i} {
	set dist [expr $s_rad2 + 10.0]
	while { $dist > [expr $s_rad2-(2.0*$sphere_rad-1.0)] } {
	    set posx [expr 2*$sphere_rad*[t_random]-$sphere_rad]
	    set posy [expr 2*$sphere_rad*[t_random]-$sphere_rad]
	    set posz [expr 2*$sphere_rad*[t_random]-$sphere_rad]
	    set dist [expr $posx*$posx + $posy*$posy + $posz*$posz]
	}
	set posx [expr $posx+[lindex $center 0] ]
	set posy [expr $posy+[lindex $center 1] ]
	set posz [expr $posz+[lindex $center 2] ] 
	part $i pos $posx $posy $posz type $type q $valency
    }
}

proc bundle_counterion_two_state_setup { n_ci sphere_rad valency {center "0 0 0"} {type "1"} {start "0"} cyl_rad cyl_length frac } {
    # Put frac*n_ci counterions inside cylinder (close to bundle)
    set n_cyl [expr int($n_ci*$frac)]
    set cyl_rad2 [expr $cyl_rad*$cyl_rad]
    for {set i $start} { $i < [expr $n_cyl+$start] } {incr i} {
	set posx [expr $cyl_length*[t_random]-($cyl_length/2.0)]
	set dist [expr $cyl_rad2 + 10.0]
	while { $dist > $cyl_rad2 } {
	    set posy [expr 2*$cyl_rad*[t_random]-$cyl_rad]
	    set posz [expr 2*$cyl_rad*[t_random]-$cyl_rad]
	    set dist [expr $posy*$posy + $posz*$posz]
	}
	set posx [expr $posx+[lindex $center 0] ]
	set posy [expr $posy+[lindex $center 1] ]
	set posz [expr $posz+[lindex $center 2] ] 
	part $i pos $posx $posy $posz type $type q $valency
    }
    # Put the rest randomly into the sphere
    set s_rad2 [expr $sphere_rad*$sphere_rad]
    for {set i [expr $n_cyl+$start]} { $i < [expr $n_ci+$start] } {incr i} {
	set dist [expr $s_rad2 + 10.0]
	while { $dist > [expr $s_rad2-(2.0*$sphere_rad-1.0)] } {
	    set posx [expr 2*$sphere_rad*[t_random]-$sphere_rad]
	    set posy [expr 2*$sphere_rad*[t_random]-$sphere_rad]
	    set posz [expr 2*$sphere_rad*[t_random]-$sphere_rad]
	    set dist [expr $posx*$posx + $posy*$posy + $posz*$posz]
	}
	set posx [expr $posx+[lindex $center 0] ]
	set posy [expr $posy+[lindex $center 1] ]
	set posz [expr $posz+[lindex $center 2] ] 
	part $i pos $posx $posy $posz type $type q $valency
    }
}

proc bundle_hair_setup { n_poly l_back cyl_rad h_dist h_length {center "0 0 0"} {type "3"} start bundle_start} {
    set mypi          3.141592653589793
    set back_id $bundle_start
    set alpha [expr -(2.0*$mypi/$n_poly)]
    for {set i 0} { $i < $n_poly } {incr i} {
	for {set j 0} { $j < $l_back } {incr j} {
	    if { $j%$h_dist == 0 } {
		set back_pos [part $back_id print pos]
		set posx [lindex $back_pos 0]
		for {set k 1} { $k <= $h_length } {incr k} {
		    set posy [expr [lindex $center 1] + ($cyl_rad-$k)*sin($i*$alpha)]
		    set posz [expr [lindex $center 2] + ($cyl_rad-$k)*cos($i*$alpha)]
		    part $start pos $posx $posy $posz type $type
		    if { $k == 1 } {
			part $start bond 0 $back_id
		    } else {
			part $start bond 0 [expr $start-1]
		    }
		    incr start
		}
	    }
	    incr back_id
	}
    }
}
