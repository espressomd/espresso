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
    }
    set usage "Usage: create_flat  \[bondl:fixz\] ]"
    array set params [::cmdline::getoptions args $options $usage]

    # First work out how many mol types there are and construct a list
    # of their lengths
    set moltypes [::mbtools::utils::listmoltypes $topology]
    set nmoltypes [llength $moltypes]
    set mollengths [::mbtools::utils::listmollengths $topology]
    set nmols [llength $topology]
    set minmoltype [::mbtools::utils::minmoltype $topology]


    set bx [lindex $boxl 0]
    set by [lindex $boxl 1]
    set bz [lindex $boxl 2]

    # We will place our lipids vertically so create a vertical
    # orientation vector
    set orient [list 0 0 1.0 ]
    

    foreach mol $topology {


	# First we choose a point in the xy plane to place the lipid the
	# standard choice is random but we might also add an option for
	# uniform placement
	lappend tailpos [expr $bx*[t_random]]
	lappend tailpos [expr $by*[t_random]]
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
	    lappend ::mbtools::system_generation::userfixedparts $i
	}
    }


     # Check particle consistency
    if { [setmd n_part] != [expr [::mbtools::utils::maxpartid $topology] + 1] } {
	mmsg::err [namespace current] "espresso has [setmd n_part] particles but [expr [::mbtools::utils::maxpartid $topology] +1] were specified in topology "
    }

}
				      


