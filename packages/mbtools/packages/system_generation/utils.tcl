#
#
# General use routines for system generation
#


namespace eval system_generation {}


#
# ::system_generation::isoutside -- 
#
# determine if a position <pos> lies outside an exclusion zone
#
# Arguments: 
# pos: The position vector to be excluded
#
# zone: A tcl list containing details about the exclusion zone.  The
# first element in the list must be the zone type. remaining elements
# will depend on the zone type.
#
proc ::system_generation::isoutside { pos zone } {
    set outside 1
    set inside 0
    switch [lindex $zone 0] {
	"sphere" {
	    set radius [lindex $zone 2]
	    set center [lindex $zone 1]
	    set dist [::mathutils::distance $center $pos]

	    if { $dist > $radius } {
		return $outside
	    } else {
		return $inside
	    }

	}
	"default" {
	    ::mmsg::err [namespace current] "[lindex $zone 0] is not a valid exclusion zone type"
	}
    }
    
}

#
# ::system_generation::maxpartid --
#  
# find the maximum particle id in a topology
# 
# 
proc ::system_generation::maxpartid { topo } { 

    set maxpart 0
    foreach mol $topo {
	for { set p 1 } { $p < [llength $mol] } { incr p } {
	    if { [ lindex $mol $p] > $maxpart } {
		set maxpart [ lindex $mol $p]
	    }
	}
    }

    return $maxpart
}

#
# ::system_generation::maxmoltypeid --
#  
# find the maximum mol type id in a topology
# 
# 
proc ::system_generation::maxmoltypeid { topo } { 
    set maxtype 0
    foreach mol $topo {
	    if { [ lindex $mol 0] > $maxtype } {
		set maxtype [ lindex $mol 0]
	    }
    }

    return $maxtype
}

#
# ::system_generation::listnmols --
#  
# Make a list containing the number of molecules of each molecule
# type.  Actually this is a list of lists. Each of the inner lists
# consists of the molecule type number and then the actual number of
# those mols
# 
# 
proc ::system_generation::listnmols { topo } { 

    set moltypes [listmoltypes $topo]
    set ntype 0

    foreach type $moltypes {
	foreach mol $topo {
	    set thistype [lindex $mol 0]
	    if { $thistype == $type } {
		incr ntype
	    }
	}       
	lappend nmolslist [list $type $ntype]
	set ntype 0
    }
    return $nmolslist
}


#
# ::system_generation::minpartid --
#
# Get the minimum particle id for this topology
# 
# 
proc ::system_generation::minpartid { topo } {  
    set startpart [lindex $topo 0 1 ]
    foreach mol $topo {
	for { set i 1 } {$i < [llength $mol] } { incr i } {
	    set thisid [lindex $mol $i]
	    if { $thisid < $startpart } {
		set startpart $thisid
	    }
	}
    }
    return $startpart
}

#
# ::system_generation::minmoltype --
#
# Get the minimum moltype for this topology
# 
# 
proc ::system_generation::minmoltype { topo } {  
    set moltypes [listmoltypes $topo]
    set startmol [lindex $moltypes 0]
    foreach tp $moltypes {
	if { $tp < $startmol } {
	    set startmol $tp
	}
    }
    return $startmol
}
#
# ::system_generation::listmoltypes --
# 
# Make a list of all the molecule types in a topology
# 
# Makes a check for duplication which would occur for an unsorted
# topology
#
proc ::system_generation::listmoltypes { topo } {  
    set currmoltype [lindex [lindex $topo 0] 0]
    foreach mol $topo {
	set moltype [lindex $mol 0]
	if { $moltype != $currmoltype } {
	    lappend typeslist $currmoltype
	    set currmoltype $moltype
	}
    }
    lappend typeslist $currmoltype

    # Check for duplication
    set reducedlist [lindex $typeslist 0]
    set duplicate 0

    for { set i 0 } { $i < [llength $typeslist] } { incr i } {
	set tp1 [lindex $typeslist $i]
	foreach tp2 $reducedlist {
	    if { $tp1 == $tp2 } {
		# we already have an entry for this type so mark it 
		set duplicate 1
	    }
	}
	if { !$duplicate } {
	    lappend reducedlist $tp1
	}
	set duplicate 0
    }


    return $typeslist
}


#
# ::system_generation::listmollengths --
# Work out the length of each molecule type and return a list of these
# lengths
#
#
proc ::system_generation::listmollengths { topo } { 
    set moltypes [listmoltypes $topo ]
 

    foreach type $moltypes {
	set lenchecksum 0
	set nchecksum 0
	foreach mol $topo {
	    set thistype [lindex $mol 0]
	    if { $thistype == $type } {
		set thislen [expr [llength $mol] - 1]
		set lenchecksum [expr $lenchecksum + $thislen]
		incr nchecksum
	    }

	}
	set lenchecksum [expr int($lenchecksum/(1.0*$nchecksum))]
	if { $lenchecksum != $thislen } {
	    mmsg::err [namespace current] "molecules of same type have different lengths"
	}
	lappend mollengths [list $type $thislen]
    }
    return $mollengths
}


#
#
#
# ::system_generation::uniquelist --
#   
# Take a list of integers and construct a list that contains no
# duplication
#
#
proc ::system_generation::uniquelist { original } { 


    lappend short [lindex $original 0]

    foreach element $original {
	if { [lsearch -exact $short $element] == -1 } {
	    lappend short $element
	}
    }

    return $short
}