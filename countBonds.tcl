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
# Last modified: 18.09.2002 by BAM                          #
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
# Created:       18.09.2002 by BAM                          #
# Last modified: 19.09.2002 by BAM                          #
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
    # total number of particles submitted to this function
    set N_T [llength $part_all]
    # Initialize the bonds-list with 'NA's
    # => if procedure is over, entries still containing 'NA' don't have any bonds
    for {set i 0} {$i<$N_T} {incr i} {
	lappend bonds "NA"
    }
    puts "Done."

    # Find the position of the bonds within 'part_all'
    puts -nonewline "        Looking for position of bond-informations... "
    for { set i 0 } {$i<$N_T} {incr i} {
	if { [findBondPos [lindex $part_all $i]]!=-1} {
	    set bond_pos [findBondPos [lindex $part_all $i]]
	    break
	}
    }
    if { $i==$N_T } {
	puts " "
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
		if { [string compare $bond_old "NA"]==0 } {
		    set bonds [lreplace $bonds $i $i "{$bond_new}"]
		} else {
		    set bonds [lreplace $bonds $i $i "$bond_old {$bond_new}"]
		}

		# ...and its respective partner.
		set bond_j [lindex [lindex $bond_i $j] 1]
		set bond_old "[lindex $bonds $bond_j]"
		set bond_new "[lindex [lindex $bond_i $j] 0] $i"
		if { [string compare $bond_old "NA"]==0 } {
		    set bonds [lreplace $bonds $bond_j $bond_j "{$bond_new}"]
		} else {
		    set bonds [lreplace $bonds $bond_j $bond_j "$bond_old {$bond_new}"]
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
