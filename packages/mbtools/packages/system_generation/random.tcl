# 
# Routines for creating a random "gas" of lipids
#

namespace eval ::system::generation {}

# ::system::generation::create_random_fluid --
#
# puts lipids in random positions with random orientation
#
proc ::system_generation::create_random_fluid {topo boxl args } {
    ::mmsg::send [namespace current] "placing lipids in a random fluid "
    set options {
	{startmol.arg     0   "starting val for molecule types"  }
	{startbdtype.arg   0  "starting val for bead types" }
	{bondl.arg     1.0   "bond length between atoms"  }
	{uniform    "use uniform lipid placement" }
	{exclude.arg "" "a region where no lipids should be placed"}
    }
    set usage "Usage: create_random_fluid topo boxl \[bondl:uniform:nhb:uniformtt:ht]"
    array set params [::cmdline::getoptions args $options $usage]
    
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

    set maxtries 1000

    
    foreach mol $topo {



	set tries 0
	set tailpos { 0 0 0 }
	set isallowed 0
	while { !$isallowed  } {

	    if {  ($tries > $maxtries) } {
		mmsg::err [namespace current] "could not place molecule exceeded max number of tries"
	    }


	    # First we choose a random point in space for the tail atom
	    lset tailpos 0 [expr $bx*[t_random]]
	    lset tailpos 1 [expr $by*[t_random]]
	    lset tailpos 2 [expr $bz*[t_random]]
	
	    if { [llength $params(exclude)] != 1 } {		
		set isallowed [isoutside $tailpos $params(exclude) ]
	    } else {
		set isallowed 1
	    }
	    incr tries

	}

	# Now choose a random orientation vector.  Actually it's not
	# strictly random if we do it like this (fix)
	lappend orient [expr [t_random]]
	lappend orient [expr [t_random]]
	lappend orient [expr [t_random]] 

	
	placemol $mol $tailpos -orient $orient -bondl $params(bondl)
	
	unset tailpos
	unset orient
    }

    return

}


























    
