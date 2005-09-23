#
#
# General use geometry type routines
#


namespace eval ::mbtools::utils {
    namespace export calc_centersofmass_bymoltype
    namespace export isoutside
    namespace export calc_com
    namespace export trap_mols
}

# ::mbtools::utils::trap_mols --
#
# Apply molecular traps to all molecules specified in a list
#
# 
#
proc ::mbtools::utils::trap_mols { molstotrap } {

    foreach mol $molstotrap {
	::mmsg::debug [namespace current] "applying trap : [lrange $mol 1 end] to mol : [lindex $mol 0]"
	analyze set trapmol [lindex $mol 0] [lindex $mol 1] [lindex $mol 2] coords [lindex $mol 3]  
    }
}

# ::mbtools::utils::calc_com --
#
# calculate center of mass of a molecule
#
# 
# Arguments:
# m : The molecule id
# topo : the topology
#
proc ::mbtools::utils::calc_com { mol } {
    set com " 0.0 0.0 0.0"
    for { set p 1 } { $p < [llength $mol ] } { incr p } { 
	set pid [lindex $mol $p]
	set pos [part $pid p p]
	for { set i 0 } { $i < 3 } { incr i } {
	    lset com $i [expr ([lindex $pos $i])/3.0 + [lindex $com $i]]
	}
    }
    return $com
}



#
# ::mbtools::utils::isoutside -- 
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
proc ::mbtools::utils::isoutside { pos zone } {
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


proc ::mbtools::utils::calc_centersofmass_bymoltype { moltypes } {
    global ::mbtools::analysis::topology
    # Extract the relevant molecules for each of the moltypes in turn
    set thismolcenters 0
    foreach moltype $moltypes {
	foreach mol $topology {	
	    if { [lindex $mol 0] == $moltype } {
		unset thismolcenters
		set cxtmp 0.0
		set cytmp 0.0
		set cztmp 0.0
		for { set i 1 } { $i < [llength $mol ] } { incr i } {
		    # Get the particle
		    set pp [part [lindex $mol $i] print pos]
		    set cxtmp [expr $cxtmp + [lindex $pp 0]]
		    set cytmp [expr $cytmp + [lindex $pp 1]]
		    set cztmp [expr $cztmp + [lindex $pp 2]]

		}
		set thisnp [expr [llength $mol] - 1]
		set thiscm [list [expr $cxtmp/($thisnp*1.0)] [expr $cytmp/($thisnp*1.0)] [expr $cztmp/($thisnp*1.0)]]
		lappend thismolcenters $thiscm	
	    }
	}
	lappend result $thismolcenters
    }
    return $result
}
