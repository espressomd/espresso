  
namespace eval ::mbtools::utils {
    namespace export maxpartid 
    namespace export maxmoltypeid
    namespace export listnmols
    namespace export minpartid
    namespace export minmoltype
    namespace export listmoltypes
    namespace export listmollengths
}

# Routines used by the analysis package which might be used by several
# different actual analysis functions


#
# ::mbtools::utils::maxpartid --
#  
# find the maximum particle id in a topology
# 
# 
proc ::mbtools::utils::maxpartid { topo } { 

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
# ::mbtools::utils::maxmoltypeid --
#  
# find the maximum mol type id in a topology
# 
# 
proc ::mbtools::utils::maxmoltypeid { topo } { 
    set maxtype 0
    foreach mol $topo {
	    if { [ lindex $mol 0] > $maxtype } {
		set maxtype [ lindex $mol 0]
	    }
    }

    return $maxtype
}

#
# ::mbtools::utils::listnmols --
#  
# Make a list containing the number of molecules of each molecule
# type.  Actually this is a list of lists. Each of the inner lists
# consists of the molecule type number and then the actual number of
# those mols
# 
# 
proc ::mbtools::utils::listnmols { topo } { 

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
# ::mbtools::utils::minpartid --
#
# Get the minimum particle id for this topology
# 
# 
proc ::mbtools::utils::minpartid { topo } {  
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
# ::mbtools::utils::minmoltype --
#
# Get the minimum moltype for this topology
# 
# 
proc ::mbtools::utils::minmoltype { topo } {  
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
# ::mbtools::utils::listmoltypes --
# 
# Make a list of all the molecule types in a topology
# 
# Makes a check for duplication which would occur for an unsorted
# topology
#
proc ::mbtools::utils::listmoltypes { topo } {  
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
# ::mbtools::utils::listmollengths --
# Work out the length of each molecule type and return a list of these
# lengths
#
#
proc ::mbtools::utils::listmollengths { topo } { 
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

