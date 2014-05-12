#############################################################
#                                                           #
# statistics.tcl                                            #
# ==============                                            #
#                                                           #
# Some scripts for statistical analysis.                    #
#                                                           #
#############################################################
#
# Copyright (C) 2010,2011,2012,2013 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#   
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
# ring_puckering
#
# Calculates the ring puckering coordinates according to
# Cremer D. and Pople J. A., J. Am. Chem. Soc. 97(6),1354-1358 (1975)
# Note that there are n-3 puckering coordinates for a n membered ring.
#
# As input you need ro give either the n coorindates of the ring atoms or
# the corresponding particle identities (internal use of part print pos).
# The parameter ret_type determines the outcome of the procedure:
# "puck" return the puckering coordinates: q_2, phi_2, ... (default) 
# "plane" return the ring plane defined by the center of mass and the unit vector perpendicular to the plane
# "disp" rerturn the n out of plane displacements z_i of the ring members
# "polar" returns spherical polar puckering coordinates (Q, theta, phi) for 6-membered rings
#
########################################################

proc ring_puckering { pos { ret_type "puck" } } {

    # number of atoms in the ring
    set n [llength $pos]

    # check for identities or positions
    if { [llength [lindex $pos 0]] == 1 } {
	set pid $pos; unset pos
	for { set i 0 } { $i < $n } { incr i } {
	    lappend pos [part [lindex $pid $i] print pos]
	}
    }

    # center of mass coordinates (here the geometrical center)
    set com { 0 0 0 }
    for { set i 0 } { $i < $n } { incr i } { 
	set com [vecadd $com [lindex $pos $i]] }
    set com [vecscale [expr 1.0/$n] $com]
    for { set i 0 } { $i < $n } { incr i } { 
	lappend pos2 [vecsub [lindex $pos $i] $com] }

    # determine ring plane
    set r1 {0 0 0}; set r2 {0 0 0}
    for { set i 0 } { $i < $n } { incr i } { 
	set r1 [vecadd $r1 [vecscale [expr sin(2*[PI]*$i/$n)] [lindex $pos2 $i] ] ]
	set r2 [vecadd $r2 [vecscale [expr cos(2*[PI]*$i/$n)] [lindex $pos2 $i] ] ]
    }
   
    set nvec [veccross_product3d $r1 $r2]
    set nvec [vecnorm $nvec]

    if { $ret_type == "plane" } { 
	set res "{$com} {$nvec}" 
	return $res
    }

    # calculate out of plane displacements
    for { set i 0 } { $i < $n } { incr i } { 
	lappend z [vecdot_product [lindex $pos2 $i] $nvec]
    }
    if { $ret_type == "disp" } { return $z }

    # calculate puckering coordinates. first pairs q_m, phi_m
    set mmax [expr ($n-1)/2]
    for { set m 2 } { $m <= $mmax } { incr m } {
	set sum1 0; set sum2 0
	for { set i 0 } { $i < $n } { incr i } { 
	    set sum1 [expr $sum1 + [lindex $z $i]*cos(2*[PI]*$m*$i/$n)]
	    set sum2 [expr $sum2 + [lindex $z $i]*sin(2*[PI]*$m*$i/$n)]
	}
	set $sum1 [expr sqrt(2.0/$n) * $sum1]
	set $sum2 [expr - sqrt(2.0/$n) * $sum2]

	set qm [expr sqrt( $sum1*$sum1 + $sum2*$sum2 ) ]
	if { $sum1 != 0 } { set phim [expr atan($sum2/$sum1)] } else { set phim [expr [PI]/2.0]}
	if { $phim < 0 } { set phim [expr [PI]+$phim] }
	lappend res $qm
	lappend res $phim
    } 
    # in case of odd number of atoms: q_(n/2)
    if { [expr $n%2] == 0 } {
	set sum 0
	for { set i 0 } { $i < $n } { incr i } { 
	    set sum [expr $sum + (pow(-1,$i)*[lindex $z $i])]
	}
	lappend res [expr $sum / sqrt($n)]
    }

    if { $ret_type == "puck" || $n != 6 } { return $res }

    # calculate spherical polar puckering coordinates for n=6 rings:
    if { $ret_type == "polar" } {
	set res2 $res; unset res
	lappend res [expr sqrt([sqr [lindex $res2 0]] + [sqr [lindex $res2 2]])]
	if { $res > 0 } { lappend res [expr acos([lindex $res2 2]/$res)] } else { lappend $res 0 }
	lappend res [lindex $res2 1]
	return $res
    }

    return "Unknown return type"
}
