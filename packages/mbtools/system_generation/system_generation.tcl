# Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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
# mbtools::system_generation --
#
# Routines in this package provide a variety of premade methods for
# setting up typical systems such a vesicles toroids and bilayers
# Author: Ira
# 

package require ::mmsg 
package require ::mbtools::utils
package provide ::mbtools::system_generation 1.0.0

# Create the empty namespace to which we shall add our routines
namespace eval ::mbtools::system_generation {
    # Global routine for setting up the system
    namespace export setup_system

    # Global variables for setup routines
    variable moltypeskey  
    variable userfixedparts
    variable trappedmols
    variable boxl
    variable middlebead
    variable topology
    variable notopo
    variable trappedmolsupdated
    variable icovermagicnums { 72 92 122 132 162 192 212 252 272 282 312 362 372 392 432 482 492 492 522 572 612 632 642 672 732 752 762 792 812 842 912 912 932 972 1002 1032 1082 1092 1112 1122 1172 1212 1242 1272 1292 1332 1332 1392 1442 1472 1472 1482 1512 1562 1572 1632 1692 1692 1712 1722 1752 1812 1832 1892 1922 1932 1962 1962 1992 2012 2082 2112 2172 2172 2192 2232 2252 2282 2292 2372 2412 2432 2442 2472 2472 2522 2562 2592 2592 2682 2712 2732 2732 2772 2792 2832 2892 2912 2922 3002 3012 3012 3042 3072 3092 3132 3162 3242 3252 3272 3312 3332 3362 3372 3432 3432 3492 3512 3612 3612 3632 3642 3642 3672 3722 3732 3792 3812 3872 3882 3972 3992 3992 4002 4032 4032 4092 4122 4172 4212 4272 4272 4322 4332 4362 4392 4412 4412 4442 4482 4532 4572 4632 4682 4692 4692 4712 4752 4812 4812 4842 4872 4892 4962 4992 5072 5072 5082 5112 5112 5132 5162 5232 5252 5292 5322 5322 5412 5432 5472 5492 5532 5532 5562 5592 5592 5672 5712 5762 5772 5792 5882 5882 5892 5892 5922 5972 6012 6032 6042 6072 6132 6192 6242 6252 6282 6312 6332 6372 6372 6372 6432 6512 6512 6522 6572 6612 6692 6732 6752 6762 6762 6792 6792 6842 6872 6882 6912 7002 7682 8192 8192 8672 9722 10832 12002 13232 14522 15872 17282 18752 20282 21872 23522 25232 27002 28832 30722 32672 34682 36752 38882 41072 43322 45632 48002 50432 52922 55472 58082 60752 63482 66272 69122 72032 75002 78032 }


    # Read in all the routines
    # Files with separate namespaces corresponding to geometries
    source [file join [file dirname [info script]] flat.tcl]
    source [file join [file dirname [info script]] random.tcl]
    source [file join [file dirname [info script]] readfile.tcl]
    source [file join [file dirname [info script]] sphere.tcl]
    source [file join [file dirname [info script]] sphere_cap.tcl]
    source [file join [file dirname [info script]] torus.tcl]
    source [file join [file dirname [info script]] cylinder.tcl]
    source [file join [file dirname [info script]] singlemol.tcl]

    # Helper function files
    source [file join [file dirname [info script]] constraint.tcl]
    source [file join [file dirname [info script]] topologies.tcl]
    source [file join [file dirname [info script]] place.tcl]

}

# ::mbtools::system_generation::setup_system -- 
#
# A large scale wrapper routine for all of the system setup commands
# in this module.
# 
# This routine should allow flexible use of multiple geometrical
# structures and topologies in a simulation.  A complete topology and
# geometry specification is provided for each object (eg sphere, flat
# bilayer etc) in the system and these are grouped into the list
# <system_specs>.  setup_system loops through each of these objects
# setting each up independently and then combining topologies at the
# end.
#
# Arguments:
# 
# system_specs: A list of specified system objects.  Each object is
# itself a list of the form < geometry n_lipidslist beads_per_mol >
# where geometry may be any allowable geometry (eg flat sphere
# etc). <n_lipidslist> is a list containing the number of molecules of
# each of the molecule types in this object and beads_per_mol
# specifies the number of particles in each of the molecule types.
#
# 
# 
#
#
proc ::mbtools::system_generation::setup_system { system_specs iboxl moltypes } {
    ::mmsg::send [namespace current] "setting up system "

    # The molecule types spec should be globally accessible
    variable moltypeskey
    variable userfixedparts
    variable trappedmols
    variable boxl
    variable topology 
    variable notopo

    set boxl $iboxl
    set moltypeskey $moltypes

    # Starting value for particle ids
    set startp 0
    # Current particle id
    set currpid $startp


    foreach spec $system_specs {
	# Flags for tracking input settings
	set geometryset 0
	set n_molslistset 0
	set n_molslist 0
	set notopo 0
	foreach item $spec {
	    switch [lindex $item 0] {
		"geometry" {
		    set geometry [lindex $item 1]
		    set geometryset 1
		}
		"n_molslist" {
		    set n_molslist [lindex $item 1]
		    set n_molslistset 1
		}
		"default" {
		    ::mmsg::warn [namespace current] "unknown item [lindex $item 0] in system spec. allowed values are: \n geometry \n n_molslist  "
		}
	    }
	}


	if { !$geometryset } {
	    mmsg::err [namespace current] "geometry not specified"
	}
	if { !$n_molslistset } {
	    mmsg::err [namespace current] "n_molslist not specified for [lindex $geometry 0]"
	}

	# Generate a topology from a list of the number and size of
	# each molecule
	foreach mtp $n_molslist {
	    set thismoltypeid [lindex $mtp 0]
	    set nmols [lindex $mtp 1]
	    set tpspec [matchtype [lindex $mtp 0]]
	    set nbeads_mol [llength [lindex $tpspec 2] ]
	    # Create the topology for this lipid type
	    set topo [create_simple_topo $nmols $nbeads_mol -moltype  $thismoltypeid -startpart $currpid ]		
		
	    # Just in case zero molecules were specified we need
	    # to check if topo was actually created at all by the
	    # last command
	    if { $topo == 0 } {
		::mmsg::err [namespace current] "no topo created for molecule type $thismoltypeid"
	    } else {
		lappend topolist $topo
		set currpid [expr [::mbtools::utils::maxpartid $topo ] + 1]
	    }

	}

	# Join all of the previously made topologies
	set first 1
	foreach topo $topolist {
	    if { $first } { 
		set topology $topo 
		set first 0
	    } else {
		set topology [join_topos $topology $topo]
	    }
	    
	}
	unset topolist


	# Now wrap the topology onto a specified geometry and perform any
	# other geometry specific tasks

	# Shuffle the topology
	# set topology [shuffle_topo $topology ]
        # NOTE: if shuffling topology is needed, use -shuffle option in each geometry

	# Now run the creation command for the specified geometry
	set createprefix "create_"
	set namespaceprefix "::mbtools::system_generation::"
	# Construct the name of the create command
	set command $geometry
	set geometry [lindex [split $geometry " "] 0]
	set createcommand "$namespaceprefix$geometry\:\:$createprefix$command "
	::mmsg::debug [namespace current] "executing $command"
	eval $createcommand
#	if { [catch  {eval $createcommand } errm ] } {
#	    mmsg::err [namespace current] "couldn't execute creation command for $command \n $errm"
#	}

	if {!$notopo} {
	    lappend topologieslist $topology
	}
	
    }

    # Join all of the previously made topologies
    set first 1
    foreach topo $topologieslist {
	if { $first } { 
	    set topology $topo 
	    set first 0
	} else {
	    set topology [join_topos $topology $topo]
	}
	
    }

    set topology [sort_topo $topology]

    return $topology

}

proc ::mbtools::system_generation::get_trappedmols {  } {
    variable trappedmols
    variable topology

    variable trappedmolsupdated

    if { [catch { set dum $trappedmols } ] } {
	::mmsg::warn [namespace current] "no trappedmols defined"
	return -1
    } else {
	if { !$trappedmolsupdated } {
	    set didntfindmol 1
	    for { set j 0 } { $j < [llength $trappedmols] } { incr j } {
		# Update trappedmols
		set fmol [lindex $trappedmols $j]
		for { set i 0 } { $i < [llength $topology] } { incr i } {
		    set mol [lindex $topology $i]
		    if { [lindex $fmol 0] == [lindex $mol 1] } {
			lset trappedmols $j 0 $i
			set didntfindmol 0
			break			
		    }
		}
		if { $didntfindmol } {
		    ::mmsg::err [namespace current] "could not get_trappedmols unable to find the corresponding particles"
		}

	    }
	}
	set trappedmolsupdated 1
	return $trappedmols
    }
}

proc ::mbtools::system_generation::get_userfixedparts {  } {
    variable userfixedparts

    if { [catch { set dum $userfixedparts } ] } {
	::mmsg::warn [namespace current] "no user fixed particles defined"
	return -1
    } else {
	return $userfixedparts
    }
}

proc ::mbtools::system_generation::get_middlebead {  } {
    variable middlebead
    return $middlebead
}
