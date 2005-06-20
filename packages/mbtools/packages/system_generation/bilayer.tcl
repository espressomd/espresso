# 
#  Routines for generating bilayer systems
#


namespace eval system_generation {}




#::system_generation::create_bilayer
#
# Place lipids in a bilayer initial configuration from a pregenerated
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
proc ::system_generation::create_bilayer {topo boxl args } {
    ::mmsg::send [namespace current] "placing lipids in a flat bilayer "

    set options {
	{startmol.arg     0   "starting val for molecule types"  }
	{startbdtype.arg   0  "starting val for bead types" }
	{uniform    "use uniform lipid placement" }
	{bondl.arg     1.0   "bond length between atoms"  }
    }
    set usage "Usage: create_bilayer topo boxl \[bondl:uniform:]"
    array set params [::cmdline::getoptions args $options $usage]
    
    variable moltypeskey

    # First work out how many mol types there are and construct a list
    # of their lengths
    set moltypes [listmoltypes $topo]
    set nmoltypes [llength $moltypes]
    set mollengths [listmollengths $topo]
    set nmols [llength $topo]
    set minmoltype [minmoltype $topo]


    set bx [lindex $boxl 0]
    set by [lindex $boxl 1]
    set bz [lindex $boxl 2]

    # We will place our lipids vertically so create a vertical
    # orientation vector
    set orient [list 0 0 1.0 ]
    

    foreach mol $topo {


	# First we choose a point in the xy plane to place the lipid the
	# standard choice is random but we might also add an option for
	# uniform placement
	lappend tailpos [expr $bx*[t_random]]
	lappend tailpos [expr $by*[t_random]]
	# We want the bilayer midplane at the box center so we need to
	# place the first tailbead half a bond length up or down
	# depending on the leaflet
	lappend tailpos [expr $bz/2.0 + 0.5*$params(bondl)*[lindex $orient 2]] 

	placemol $mol $tailpos -orient $orient -bondl $params(bondl)

	#  Since we are just making a bilayer we alternate between upper
	#  and lower layers and therefore orient should be reversed
	set orient [::mathutils::scalevec $orient -1]
	
	unset tailpos
	
    }
    

}
				      


