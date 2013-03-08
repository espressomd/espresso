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
# Proceedures for analyzing the lipid flip-flop rate
#

# ::mbtools::analysis::Compliment --
#
# Function called by "analyze_flipflop".  Compliment takes two lists
# of particle orientations and compares them to find the number of
# particles whose orientation is different between the two.
#

namespace eval ::mbtools::analysis {}

namespace eval ::mbtools::analysis::flipflop {
    variable av_flip 0.0
    variable av_flip_i 0
    variable f_tvsflip
    variable l_orients
    variable l_orients_start
    variable verbose

    namespace export setup_flipflop
    namespace export analyze_flipflop
    namespace export printav_flipflop
    namespace export resetav_flipflop

}

proc ::mbtools::analysis::flipflop::resetav_flipflop { } {
    variable av_flip 
    variable av_flip_i 
    set av_flip 0
    set av_flip_i 0
}


proc Compliment { list1 list2 } {
    set flip 0
    set bilcnt 0

    for { set i 0 } { $i < [llength [lindex $list1 1] ] } {incr i } { 
	#	    puts "[lindex $list1 1 $i ] : [lindex $list2 1 $i]"
	if { [lindex $list1 1 $i ] != [lindex $list2 1 $i] && [lindex $list1 1 $i ] != 2 && [lindex $list2 1 $i] != 2 && [lindex $list1 1 $i ] != 3 && [lindex $list2 1 $i] != 3 } {
	    incr flip
	    
	    #	puts "flip $flip : [lindex $list1 1 $i ] : [lindex $list2 1 $i] "
	    #	flush stdout
	}
	if { [lindex $list1 1 $i ] != 2 && [lindex $list2 1 $i] != 2 && [lindex $list1 1 $i ] != 3 && [lindex $list2 1 $i] != 3 } {
	    incr bilcnt
	}
    }
    if { $bilcnt == 0 } {
	::mmsg::warn [namespace current] "no lipids in bilayer: cannot perform flip-flip analysis"
	puts "$list1 list2"
	return 0
    }
    return [expr ($bilcnt - $flip)/(1.0*$bilcnt)]
    
}

proc ::mbtools::analysis::flipflop::printav_flipflop { } {
    variable f_tvsflip
    variable av_flip
    variable av_flip_i
    global ::mbtools::analysis::time
    if { $av_flip_i > 0 } {
	puts $f_tvsflip "$time [expr $av_flip/(1.0*$av_flip_i)]"
    }
    flush $f_tvsflip
}

proc ::mbtools::analysis::flipflop::setup_flipflop { args } {
    global ::mbtools::analysis::outputdir
    global ::mbtools::analysis::suffix
    global ::mbtools::analysis::iotype
    variable f_tvsflip
    variable verbose
    variable l_orients_start
    set options {
	{verbose "print out lots of stuff" }
    }
    set usage "Usage: setup_flipflop gridm:straycutoff:verbose "
    array set params [::cmdline::getoptions args $options $usage]

    set verbose $params(verbose)

    ::mmsg::debug [namespace current] "opening $outputdir/time_vs_flip$suffix "

		if { [file exists "$outputdir/time_vs_flip$suffix"] } {
		    set newfile 0
		} else { 
		    set newfile 1
		}

		set f_tvsflip [open "$outputdir/time_vs_flip$suffix" $iotype]
		set f_tvsflip [open "$outputdir/time_vs_flip$suffix" $iotype]
		if { $newfile || $iotype == "w"} {
		    puts $f_tvsflip "\# Note: F(t) = N_up(t) + N_d(t)/N. See long membrane paper for details"
		    puts $f_tvsflip "\# Time flip_parameter"
		}
		set l_orients_start [ analyze get_lipid_orients ]

}


# ::mbtools::analysis::analyze_flipflop --
#
#  Obtains the flip state of the current system and compares it to the
#  starting state to determine the number of flips that have occurred.
#
proc ::mbtools::analysis::flipflop::analyze_flipflop { } {
    variable av_flip
    variable av_flip_i
    variable f_tvsflip
    variable l_orients
    variable l_orients_start
    variable verbose

    ::mmsg::send [namespace current] "analyzing flip-flop"

    set l_orients [ analyze get_lipid_orients ]
    
    set flip [ Compliment $l_orients $l_orients_start ]
    set av_flip [expr $av_flip + $flip ]
    incr av_flip_i
    
    if { $verbose } { 
	::mmsg::send [namespace current] "flipflip: $flip <flipflop>: [expr $av_flip/($av_flip_i*1.0) ]"
	flush stdout
    }
}


