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
#  Routines for generating flat systems
#


namespace eval mbtools::system_generation {}



namespace eval ::mbtools::system_generation::flat {


    namespace export create_flat
}

#::mbtools::system_generation::flat::create_flat
#
# Place lipids in a flat initial configuration from a pregenerated
# topology.
#
# This routine is really only designed to work for simple linear
# lipids.  It supports topologies with unlimited lipid types provided
# that all have the same number of head beads. 
#
#
# Arguments: 
#
# topo: The topology in espresso format. It should be
# sorted in ascending order according to molecule type (not bead id).
#
# Options:
# ht: All head beads are type 0 and all tails are uniform in type
# uniformtt: All tails are created from just one bead type (not yet
# implemented
#
#
#
proc ::mbtools::system_generation::flat::create_flat {args } {
    ::mmsg::send [namespace current] "placing lipids in a flat bilayer"

    global ::mbtools::system_generation::moltypeskey
    global ::mbtools::system_generation::boxl
    global ::mbtools::system_generation::topology

    set options {
	{bondl.arg     1.0   "bond length between atoms"  }
	{fixz "fix the z positions of all particles in this bilayer"}
	{crystal "set the lipids on a grid, instead of a random placement"}
	{half "fills only half ot the xy plane with lipids, PBC in only one direction"}
	{pancake "create a pancake bilayer, does not extend through the PBC"}
        {shuffle "shuffle topology before placement"}
    }
    set usage "Usage: create_flat  \[bondl:fixz:crystal:half:pancake:shuffle\] ]"
    array set params [::cmdline::getoptions args $options $usage]

    # First work out how many mol types there are and construct a list
    # of their lengths
    set moltypes [::mbtools::utils::listmoltypes $topology]
    set nmoltypes [llength $moltypes]
    set mollengths [::mbtools::utils::listmollengths $topology]
    set nmols [llength $topology]
    set minmoltype [::mbtools::utils::minmoltype $topology]

    if { $params(shuffle) } {
        set topology [::mbtools::system_generation::shuffle_topo $topology ]
    }

    set bx [lindex $boxl 0]
    set by [lindex $boxl 1]
    set bz [lindex $boxl 2]

    # We will place our lipids vertically so create a vertical
    # orientation vector
    set orient [list 0 0 1.0 ]
    
    # need dummy variables to iterate for the crystal structure
    set i 0; set j 0
    # ratio of two lengths r
    set r [expr $bx/$by]
    # number of molecules
    set Npart 0
    foreach mol $topology {
	incr Npart
    }
    # determine the number of lipids along each axis
    set N_x [expr sqrt($Npart*$r)]
    set N_y [expr sqrt($Npart/$r)]

    if ($params(pancake)) {
	# assume area per lipid of 1.25, then determine radius from the number
	# of lipids on ONE leaflet. Return an error if the area of the pancake
	# is larger than the area of the box.
	set area_per_lipid 1.25
	set d_bin [expr { sqrt(1.25) }]
	set radius_pancake [expr sqrt($area_per_lipid*$Npart/2./3.14159265)]
	set n_bin [expr { ceil($radius_pancake/$d_bin) }]
	if {pow($radius_pancake,2)*3.14159265 > $bx*$by} {
	    mmsg::err [namespace current] "Area of pancake bilayer is larger than the box size"
	}
	# find the center of the box
	set center_x [expr {0.5 * $bx}]
	set center_y [expr {0.5 * $by}]
    }

    foreach mol $topology {
	# CRYSTAL STRUCTURE
	if ($params(crystal)) {
	    incr i
	    if {$i*$bx/$N_x>$bx} { 
		set i 0
		incr j
	    }
	    lappend tailpos [expr $i*$bx/$N_x]
	    lappend tailpos [expr $j*$by/$N_y]
	} elseif ($params(half)) { 
	    # Half bilayer
	    lappend tailpos [expr .25*$bx+.5*$bx*[t_random]]
	    lappend tailpos [expr $by*[t_random]]
	} elseif {$params(pancake)} {
	    # Create a bounded pancake (does not extend through the PBC)
	    # in the x-y plane. Make sure the pancake is smaller than the box
	    if {$radius_pancake > $bx || $radius_pancake > $by} {
		mmsg::err [namespace current] "trying to create a pancake bilayer that's larger than the box size"
	    }
	    # Position all lipids randomly
	    # set circle_r [expr $radius_pancake * [t_random]]
	    # set circle_theta [expr 6.283185307* [t_random]];# 2*pi
	    # lappend tailpos [expr $circle_r * cos($circle_theta) + $bx/2.]
	    # lappend tailpos [expr $circle_r * sin($circle_theta) + $by/2.]
	   
	    # Position all lipids on lattice
	    incr i
	    if {$i>[expr { $n_bin*2 } ] } {
		set i 0
		incr j
	    }
	    # Since it is a square lattice we need to check if the node is indeed in our pancake
	    set center_distance [expr sqrt( ($n_bin - $i) * ($n_bin - $i) + ($n_bin - $j) * ($n_bin - $j)) * $d_bin ]	    
	    while { $center_distance > $radius_pancake } { 
		incr i
		if {$i>[expr { $n_bin*2 } ] } {
		    set i 0
		    incr j
		}
		if {$j>[expr { $n_bin*2 } ] } { set j 0}
		set center_distance [expr sqrt( ($n_bin - $i) * ($n_bin - $i) + ($n_bin - $j) * ($n_bin - $j)) * $d_bin ]
	    }
	    set tmp_x [expr $center_x - $d_bin * $n_bin + $i * $d_bin]
	    set tmp_y [expr $center_y - $d_bin * $n_bin + $j * $d_bin]
	    lappend tailpos $tmp_x
	    lappend tailpos $tmp_y
	   
 	} else {
 	    # RANDOM PLACEMENT
 	    # First we choose a point in the xy plane to place the lipid the
 	    # standard choice is random but we might also add an option for
 	    # uniform placement
 	    lappend tailpos [expr $bx*[t_random]]
 	    lappend tailpos [expr $by*[t_random]]
 	}
 	# We want the flat midplane at the box center so we need to
 	# place the first tailbead half a bond length up or down
 	# depending on the leaflet
 	lappend tailpos [expr $bz/2.0 + 0.5*$params(bondl)*[lindex $orient 2]] 
	
 	::mbtools::system_generation::placemol $mol $tailpos -orient $orient -bondl $params(bondl)
	
 	#  Since we are just making a flat we alternate between upper
 	#  and lower layers and therefore orient should be reversed
 	set orient [::mbtools::utils::scalevec $orient -1]
		
	unset tailpos
    }
    
    if ($params(fixz)) {
	for {set i [::mbtools::utils::minpartid $topology] } { $i <  [setmd n_part] } {incr i} {
	    set fixvalue [part $i print fix]
	    #check to see if the particle has been otherwise specially fixed
	    part [expr $i] fix [lindex $fixvalue 0] [lindex $fixvalue 1] 1
	    #	    lappend ::mbtools::system_generation::userfixedparts $i
	}
    }


    # Check particle consistency
    if { [setmd n_part] != [expr [::mbtools::utils::maxpartid $topology] + 1] } {
	mmsg::err [namespace current] "espresso has [setmd n_part] particles but [expr [::mbtools::utils::maxpartid $topology] +1] were specified in topology "
    }
    
}



