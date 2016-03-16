# Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#   Max-Planck-Institute for Polymer Research, Theory Group
#  
# This file is part of ESPResSo.
#  
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 
#  
#############################################################
#                                                           #
# countBonds.tcl                                            #
# ==============                                            #
#                                                           #
# countBonds                                                #
# ----------                                                #
# Analyzes the topology of a given configuration of         #
# particles, returning a tcl-list containing all bonds      #
# leading to or emerging from the respective particle.      #
#                                                           #
# Input:                                                    #
# A tcl-list with informations about each particle,         #
# formated the same way as the output of '[part]'           #
#                                                           #
# Output:                                                   #
# A tcl-list containing the bonds for each particle,        #
# formated similarly to the bond-section in the output of   #
# '[part]', i.e. '{ {type partner} {type partner} ... }'    #
#                                                           #
# Created:       17.09.2002 by BAM                          #
# Last modified: 26.09.2002 by BAM                          #
#                                                           #
#                                                           #
# findBondPos                                               #
# -----------                                               #
# Analyzes the parameter set of a given particle to find    #
# out at which position the informations about the bonds    #
# are stored, returning that index.                         #
#                                                           #
# Input:                                                    #
# A tcl-list of parameters of a particle, formated the      #
# same way as the output of '[part <part-num>]'             #
#                                                           #
# Output:                                                   #
# The index of the bond-list such that                      #
# '[lindex <input> <output>]' will return the list of bonds #
# '{<type> <partner>} {<type> <partner>} ...';              #
# returns '-1' if no bonds have been found in <input>.      #
#                                                           #
# Created:       18.09.2002 by BAM                          #
# Last modified: 18.09.2002 by BAM                          #
#                                                           #
#                                                           #
# findPropPos                                               #
# -----------                                               #
# Analyzes the parameter set of a given particle to find    #
# out at which position the informations about a certain    #
# property is stored.                                       #
# This function generalizes 'findBondPos', i.e. using it    #
# with 'bonds' as 2nd argument returns the same results.    #
#                                                           #
# Input:                                                    #
# - A tcl-list of parameters of a particle, formated the    #
#   same way as the output of '[part <part-num>]'           #
# - The name of the property for which the position should  #
#   be found                                                #
# e.g. 'findPropPos [part 13] bonds' (usually returns '17') #
#                                                           #
# Output:                                                   #
# The index where the properties are stored,                #
# e.g. if 'bonds' are the 2nd argument the index of the     #
# bonds is returned such that '[lindex <input> <output>]'   #
# will yield the list of bonds                              #
# '{<type> <partner>} {<type> <partner>} ...'.              #
# However, '-1' is returnd if <property> is not in <input>. #
#                                                           #
#                                                           #
#############################################################
proc findBondPos { tmp_part } {
    for {set i 0} {$i<[llength $tmp_part]} {incr i} {
	if { [string compare [lindex $tmp_part $i] "bonds"]==0 } {
	    return [expr $i+1]
	}
    }
    return -1
}   




proc findPropPos { tmp_part tmp_prop} {
    for {set i 0} {$i<[llength $tmp_part]} {incr i} {
	if { [string compare [lindex $tmp_part $i] $tmp_prop]==0 } {
	    return [expr $i+1]
	}
    }
    return -1
}   




proc countBonds { part_all } {
    puts "    Analysis of topology in progress..."

    puts -nonewline "        Preparing environement... "
    flush stdout
    # total number of particles submitted to this function
    set N_T [llength $part_all]
    # Initialize the bonds-list with the particle_number (taken from $part_all)
    # => if procedure is over, entries of llength==1 don't have any bonds
    for {set i 0} {$i<$N_T} {incr i} {
	lappend bonds [lindex [lindex $part_all $i] 0]
    }
    puts "Done."

    # Find the position of the bonds within 'part_all'
    puts -nonewline "        Looking for position of bond-informations... "
    flush stdout
    for { set i 0 } {$i<$N_T} {incr i} {
	if { [findBondPos [lindex $part_all $i]]!=-1} {
	    set bond_pos [findBondPos [lindex $part_all $i]]
	    break
	}
    }
    if { $i==$N_T } {
	puts "Failed.\n"
	puts "----------------"
	puts "--> Warning! <--"
	puts "----------------"
	puts "--> None of the $N_T particles have any bonds!"
	puts "Aborting...\n"
	exit
    }
    puts "found them at index $bond_pos in particle $i."
    

    puts -nonewline "        Analyzing and adding bonds now... "
    flush stdout
    set bond_nr 0
    set old [expr 10*$N_T]
    for {set i 0} {$i<$N_T} {incr i} {
	# Draw a status bar for every 10% of progress
	if { [expr $i*100]>=$old } {
	    set old [expr $old + 10*$N_T]
	    puts -nonewline ".,."
	    flush stdout
	}

	# Informations about the bonds is usually stored in slots '[expr $bond_pos-1]' & '$bond_pos'
	# of the particle with the lower index number, but not at the bonding partner.
	# Hence, these informations are transferred to 'bonds' at both the particle's and its partner's index.
	if { [llength [lindex $part_all $i]>[expr $bond_pos -1]] } {
	    set bond_i [lindex [lindex $part_all $i] $bond_pos]
	    set bond_l [llength [lindex [lindex $part_all $i] $bond_pos]]

	    # Look at all '$bond_l' partners this $i-th particle has 
	    # (stored in $bond_i containing all the bonding informations about $i)
	    for {set j 0} {$j<$bond_l} {incr j} {

		# Add bonding informations to the particle...
		set bond_old "[lindex $bonds $i]"
		set bond_new "[lindex $bond_i $j]"
		set bonds [lreplace $bonds $i $i "$bond_old {$bond_new}"]

		# ...and its respective partner.
		set bond_j [lindex [lindex $bond_i $j] 1]
		set bond_j_f -1
		for {set k 0} {$k<$N_T} {incr k} {
		    # since the particle_number does not have to correspond to the index of that particle in $part_all
		    # (e.g. particle 45 might be stored at [lindex $part_all 113], for whatever strange reason)
		    # one has to find the tcl-list-index for $bond_j to know where to save its bonding informations
		    if { [lindex [lindex $part_all $k] 0]==$bond_j } {
			set bond_j_f $k
			break
		    }
		}
		if { $bond_j_f > -1 } {
		    set bond_j $bond_j_f
		    set bond_old "[lindex $bonds $bond_j]"
		    set bond_new "[lindex [lindex $bond_i $j] 0] [lindex [lindex $part_all $i] 0]"
		    set bonds [lreplace $bonds $bond_j $bond_j "$bond_old {$bond_new}"]
		} else {
		    puts "Failed.\n"
		    puts "----------------"
		    puts "--> Warning! <--"
		    puts "----------------"
		    puts "--> Cannot find an entry for particle $bond_j in the supplied list of $N_T particles!"
		    puts -nonewline "--> (Got:"
		    for {set k 0} {$k<$N_T} {incr k} { puts -nonewline " [lindex [lindex $part_all $k]]" }
		    puts ")"
		    puts "Aborting...\n"
		    exit
		}

		# Count the amount of bonds set.
		incr bond_nr
	    }
	}
    }
    puts ".,. got $bond_nr bonds total."


    # Now return what has been found
    puts "    Analysis completed. Returning control to calling script..."
    return $bonds
}
