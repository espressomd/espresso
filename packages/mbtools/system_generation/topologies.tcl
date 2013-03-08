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
#  Proceedures for creating various types of topologies.  In this case
#  topology refers to the espresso feature with the same name and
#  describes the groupings of atoms into molecules.
#


namespace eval mbtools::system_generation {}


#----------------------------------------------------------#
# ::mbtools::system_generation::create_simple_topo --
#
# Creates a basic topology from a single molecule type
#
proc ::mbtools::system_generation::create_simple_topo {n_mols beads_per_mol  args } {

    set options {
	{moltype.arg 0 "type of molecule to create" }
	{startpart.arg 0 "The starting id for particles" }
    }
    set usage "Usage: create_simple_topo n_mols beads_per_mol \[startmol:startpart:]"
    array set params [::cmdline::getoptions args $options $usage]

    ::mmsg::send [namespace current] "creating a topology for moltype $params(moltype)" 

    # First we create a very simple completely ordered and sorted
    # topology 
    set mol_type $params(moltype)
    set partid $params(startpart)

    if { $n_mols <= 0 } {
	return 0
    }

    for { set i 0 } { $i < $n_mols } { incr i } {
	# construct the topo_entry for this lipid

	set topo_entry $mol_type
	for { set j 0 } { $j < $beads_per_mol } { incr j } {
	    lappend topo_entry $partid
	    incr partid
	}
	lappend topo $topo_entry
    }
    return $topo

}



# 
# ::mbtools::system_generation::shuffle_topo--
#  
# Performs a shuffle on the topology so that molecules are no longer
# in blocks of the same type
#
#
# Arguments:
#
# topo: The topology to be shuffled.  
# 

proc ::mbtools::system_generation::shuffle_topo {topo} {
    mmsg::send [namespace current] "shuffling topology" 

    set startmol [::mbtools::utils::minmoltype $topo ]
    set startpart [::mbtools::utils::minpartid $topo ]

    set full_n_molslist [::mbtools::utils::listnmols $topo]
    set n_molstotal [llength $topo]

    set remaining_molslist $full_n_molslist

    # Check that we have there are no discontinuities in the particle ids

    # To do the shuffle properly we really need to completely regenerate the topology


    # The first step is to construct a list with randomly allocated
    # molecule types to serve as a template
    set n_remaining $n_molstotal
    mmsg::send [namespace current] "constructing topo template " nonewline 
    flush stdout
    set dotfreq [expr int(floor($n_molstotal/10.0))]
    if { $dotfreq < 1 } { set dotfreq 1}
    for { set i 0 } { $i < $n_molstotal } { incr i } {
	
	# According to the molecule proportions determine the type of
	# the next molecule
	set lpick [expr $n_remaining*[t_random]]
	for {set lnum 0 } { $lnum <  [llength $remaining_molslist] } { incr lnum } {
	    set lpick [expr $lpick - [lindex $remaining_molslist $lnum 1]]

		if  { $lpick <= 0 } {
		    set thislipidtype [lindex $remaining_molslist $lnum 0]
		    # Subtract the picked lipid from our molslist
		    set n_before [lindex  $remaining_molslist $lnum 1]
		    lset remaining_molslist $lnum 1 [expr  $n_before -1]
		    set n_remaining [expr $n_remaining -1 ]
		    break
		}
	}

	lappend typetemplate $thislipidtype

        if { $i%$dotfreq ==0 } {
	    mmsg::send [namespace current] "." nonewline
	    flush stdout
	}
    }
    mmsg::send [namespace current] "done" 

    # We need a list telling us how many atoms are in each molecule
    set molsizes [::mbtools::utils::listmollengths $topo ]
    set pnum $startpart
    foreach type $typetemplate {
	foreach t $molsizes {
	    if { [lindex $t 0] == $type  } {
		set thislen [lindex $t 1]
	    }
	}
	set thismol $type
	for { set p 0 } { $p < $thislen } { incr p } {
	    lappend thismol $pnum
	    incr pnum
	}
	lappend shuffledtopo $thismol
    }
 



    flush stdout
    return $shuffledtopo
}

#
# ::mbtools::system_generation::sort_topo -- 
#
# Sort a topology into ascending molecule type order
#
#
proc ::mbtools::system_generation::sort_topo { topo } {    
    set n_molstotal [llength $topo]
    set maxtp [::mbtools::utils::maxmoltypeid $topo]

    for { set l 0 } { $l <= $maxtp } { incr l } {
	for { set i 0 } { $i < $n_molstotal } { incr i } {
	    if { [lindex [lindex $topo $i] 0 ] == $l } {
		lappend sortedtopo [lindex $topo $i]
	    }
	}
    }
    return $sortedtopo
}


#
# ::mbtools::system_generation::join_topos--
#
# Join two topologies
#
proc ::mbtools::system_generation::join_topos { topo1 topo2 args } {
    mmsg::send [namespace current] "joining topologies"

    set options {
	{sortmolid  "sort by molid" }
    }
    set usage "Usage: join_topos \[sortmolid]"
    array set params [::cmdline::getoptions args $options $usage]
    foreach mol $topo1 {
	lappend joinedtopo $mol
    }
    foreach mol $topo2 {
	lappend joinedtopo $mol
    }

    return $joinedtopo
}