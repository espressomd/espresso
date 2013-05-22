# Copyright (C) 2010,2012,2013 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#   Max-Planck-Institute for Polymer Research, Theory Group
#  
# This file is part of ESPResSo.
#   
# ESPResSo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#  
# ESPResSo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
# This set of routines constructs a spherical cap where we attempt
# to place the lipids at a uniform density taking account of the
# sphere's curvature. 
#

namespace eval ::mbtools::system_generation {}

namespace eval ::mbtools::system_generation::sphere_cap {


    namespace export create_sphere_cap
}

#
# ::mbtools::system_generation::sphere::create_sphere_cap --
#
# Place a topology into a spherical cap geometry. This cap is part of a sphere
# with radius r as an option. The area of the cap depends on the number of molecules
# in the topology and the option "initarea". This is done by placing molecules
# from the polar till all molecules in topology have been placed.
# 
# options:
#         half  Create a half of a sphere, whose radius is calculated by
#               the number of molecules and initarea.

proc ::mbtools::system_generation::sphere_cap::create_sphere_cap { args } { 
    
    # ---- Process Command Line Args --------------------------- #
    
    set options {
	{r.arg   10.0       "the radius of the sphere" }
	{half  "create a half sphere" }
	{c.arg  { 0.0 0.0 0.0 } "location of the center of the sphere " }
	{initarea.arg      1.29    "the starting value for area per mol" }
	{shuffle "shuffle the topology before placing molecules "}
	{bondl.arg     1.0   "bond length between atoms"  }
    }
    set usage "Usage: create_sphere_cap \[r:half:c:initarea:shuffle:bondl\] ]: "
    array set params [::cmdline::getoptions args $options $usage]
    
    global ::mbtools::system_generation::moltypeskey
    global ::mbtools::system_generation::boxl
    global ::mbtools::system_generation::topology

    puts "$params(shuffle)"   
  
    if { $params(shuffle) } {
	set topology [::mbtools::system_generation::shuffle_topo $topology ]
    }


    #Construct the center
    set center [list 0.0 0.0 0.0 ]
    lset center 0 [expr [lindex $params(c) 0]]
    lset center 1 [expr [lindex $params(c) 1]]
    lset center 2 [expr [lindex $params(c) 2]]

    # First work out how many mol types there are and construct a list
    # of their lengths

    set moltypes [::mbtools::utils::listmoltypes $topology]
    set nmoltypes [llength $moltypes]
    set mollengths [::mbtools::utils::listmollengths $topology]
    set nmols [llength $topology]
    set minmoltype [::mbtools::utils::minmoltype $topology]
   
    set pancake_area [expr 0.5*$nmols*$params(initarea)]
    # calculate the radius of our sphere
    if ($params(half)) {
	set rad_av [expr sqrt($pancake_area/(2*3.14159))]
	mmsg::send [namespace current] "Creating a half sphere."
    } else {
	set rad_av [expr $params(r)*1.0]
    }
    mmsg::send [namespace current] "placing lipids on a spherical cap of radius $rad_av "

    # The surface area of the sphere should be larger than or equal to the area of pancake, thus a lower
    # bound for r exists.
    set r_lowerbound [expr sqrt($nmols*$params(initarea)/(8.0*3.14159))]
    if { $params(r) < $r_lowerbound } {
	mmsg::err [namespace current] "Radius is too small to hold $nmols molecules. Input a number larger than $r_lowerbound."
    }
    # Work out the average molecule size
    set sum 0
    foreach mol $topology {
	set sum [expr $sum + [llength $mol] - 1]
    }
    set avmolsize [expr $sum/(1.0*$nmols)]

    # Make some aliases for useful numbers
    set del [expr $avmolsize/2.0]
    set pi 3.141592
    set marea $params(initarea)

    

    #Based on the given radius and average molecule size work out
    #how many should go on outer and inner leaflets
	set ratio_r [expr ($rad_av-$del)/($rad_av+$del)]
	set ratio_n [expr $ratio_r*$ratio_r]
    set n_mols_outer [expr $nmols*1.0/(1.0+$ratio_n)]
    set n_mols_inner [expr $n_mols_outer*$ratio_n]
    mmsg::send [namespace current] "n+: $n_mols_outer n-:$n_mols_inner"
    set n_mols_outer [expr int(($n_mols_outer))]
    set n_mols_inner [expr $nmols - $n_mols_outer]
    mmsg::send [namespace current] "n+: $n_mols_outer n-:$n_mols_inner"


    # Determine radius for the outer layer
    set rad [expr $rad_av + $del ]
    mmsg::send [namespace current] "outer radius: $rad"


    # --- Outer Layer -------- #
    # Now go through and place all mols in the outer layer
    set d [expr sqrt($marea)]
    set mtheta [expr int(($pi*$rad/$d))]
    set dtheta [expr ($rad*$pi)/(1.0*$mtheta)]
    set dphi [expr $marea/$dtheta]

    set molnum 0
    set outercount 0
    set removecount 0
    for { set m 0 } { $m < $mtheta } { incr m } {
	set theta [expr $pi*($m+0.5)/($mtheta*1.0)]
	set mphi [expr int((2.0*$pi*$rad*sin($theta)/$dphi))]
	for { set n 0 } { $n < $mphi } { incr n } {
	    set phi [expr 2*$pi*$n/(1.0*$mphi)]
	    # Find the position of the first tail bead in the molecule
	    lappend tailpos [expr ($rad_av+0.5)*sin($theta)*cos($phi) + [lindex $center 0]]
	    lappend tailpos [expr ($rad_av+0.5)*sin($theta)*sin($phi) + [lindex $center 1]]
	    lappend tailpos [expr ($rad_av+0.5)*cos($theta) + [lindex $center 2] ]
	    
	    # Find the orientation vector for the molecule that is
	    # normal to the sphere
	    lappend orient [expr sin($theta)*cos($phi)]
	    lappend orient [expr sin($theta)*sin($phi)]
	    lappend orient [expr cos($theta)]


	    # Place the molecule
	    set mol [lindex $topology $molnum]
	    ::mbtools::system_generation::placemol $mol $tailpos -orient $orient -bondl $params(bondl)

	    incr molnum
	    incr outercount
	   
	    unset tailpos
	    unset orient
	    if { $outercount == $n_mols_outer } {
		break
	    }
	}

	if { $outercount == $n_mols_outer } {
		break
	}
    }
    mmsg::send [namespace current] "placed [expr $outercount] mols in outer layer"
   #mmsg::send [namespace current] "$removecount mols were skipped to get this number"
    
    #------------------ Inner ---------#
    # Determine r1 and r2 for the inner layer
    set rad [expr $rad_av - $del ]
    mmsg::send [namespace current] "inner radius: $rad "
    
    
    # Now go through and place all mols in the Inner layer
    set d [expr sqrt($marea)]
    set mtheta [expr int(($pi*$rad/$d))]
    set dtheta [expr ($rad*$pi)/(1.0*$mtheta)]
    set dphi [expr $marea/$dtheta]
    set innercount 0
    set removecount 0
    for { set m 0 } { $m < $mtheta } { incr m } {
	set theta [expr $pi*($m+0.5)/($mtheta*1.0)]
	set mphi [expr int((2.0*$pi*$rad*sin($theta)/$dphi))]
	for { set n 0 } { $n < $mphi } { incr n } {
	    set phi [expr 2*$pi*$n/(1.0*$mphi)]
	    lappend tailpos [expr ($rad_av-0.5)*sin($theta)*cos($phi) + [lindex $center 0]]
	    lappend tailpos [expr ($rad_av-0.5)*sin($theta)*sin($phi) + [lindex $center 1]]
	    lappend tailpos [expr ($rad_av-0.5)*cos($theta) + [lindex $center 2] ]

	    # Find the normal vector pointing inwards
	    lappend orient [expr -sin($theta)*cos($phi)]
	    lappend orient [expr -sin($theta)*sin($phi)]
	    lappend orient [expr -cos($theta)]
   
    
	    # Place the molecule
	    set mol [lindex $topology $molnum]
	    ::mbtools::system_generation::placemol $mol $tailpos -orient $orient -bondl $params(bondl)
		
	    incr molnum
	    incr innercount
	    
	    unset tailpos
	    unset orient
	    if { $innercount == $n_mols_inner } {
		break
	    }
	}
	if { $innercount == $n_mols_inner } {
	    break
	}
    }


    # Check particle consistency
    if { [setmd n_part] != [expr [::mbtools::utils::maxpartid $topology] + 1] } {
	mmsg::err [namespace current] "espresso has [setmd n_part] particles but [expr [::mbtools::utils::maxpartid $topology] +1] were specified in topology "
    }

    mmsg::send [namespace current] "placed [expr $innercount] mols in inner layer"
    #mmsg::send [namespace current] "$removecount random mols were skipped to get this number"
    mmsg::send [namespace current] "total mols in vesicle: $molnum"
    mmsg::send [namespace current] "uniform spherical vesicle created" 
    flush stdout

    return
}

