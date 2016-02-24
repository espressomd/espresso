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
#
# Procedures for determining the distribution of local lipid orientations
#

namespace eval ::mbtools::analysis {}


namespace eval ::mbtools::analysis::localorients {
    variable av_localorients 0
    variable av_localorients_i 0
    variable verbose
    variable nbins
    variable range

    namespace export setup_localorients
    namespace export analyze_localorients
    namespace export printav_localorients
    namespace export resetav_localorients
}

proc ::mbtools::analysis::localorients::resetav_localorients { } {
    # Do nothing
}

proc ::mbtools::analysis::localorients::printav_localorients {} {
    variable av_localorients_i
    variable av_localorients
    variable f_localorients
    global ::mbtools::analysis::outputdir
    global ::mbtools::analysis::suffix
    variable nbins
    variable range

    if { $av_localorients_i > 0 } {
		    
	set f_localorients [open "$outputdir/av_localo$suffix" "w" ]
	# Write a header to the file
	puts $f_localorients "\# local orientation distribution of lipids "
	puts $f_localorients ""

		   
	set lobinwidth [expr $range/($nbins*1.0)]

	# Write out the local orientation distribution to a file
	for { set bin 0 } { $bin < [llength $av_localorients] } { incr bin } {
	    set currbin [lindex $av_localorients $bin]
	    puts $f_localorients "[expr $bin*$lobinwidth] [expr $currbin/(1.0*[llength $topology])]"
	}		    
	close $f_localorients
    }

}


proc ::mbtools::analysis::localorients::setup_localorients { args } {
    variable verbose
    variable range
    variable nbins
    variable av_localorients

    set options {
	{verbose "print out lots of stuff" }
	{range.arg "1.0" "Range for localorients calculation" }
	{nbins.arg "100" "Number of bins used for analysis" }
    }
    set usage "Usage: setup_pressure verbose:range "
    array set params [::cmdline::getoptions args $options $usage]
    
    set verbose $params(verbose)
    set range $params(range)
    set nbins $params(nbins)

    set av_localorients 0
    unset av_localorients
    #Initialize av_localorients
    for { set bn 0 } { $bn < $nbins } { incr bn } {
	lappend av_localorients 0
    }
}


# ::mbtools::analysis::analyze_localorients --
#  
# Calculates the projection of the lipid orientation vector onto the
# xy plane for each lipid and then bins the absolute values of these
# vectors.
#
proc ::mbtools::analysis::localorients::analyze_localorients { printflag } {

    variable av_localorients
    variable av_localorients_i
    variable nbins
    variable range
    ::mmsg::send [namespace current] "analyzing local orients"

    set nbins $nbins
    set binwidth [expr $range/($nbins*1.0)]

    # Internal copy of l_orients ( same molecule order as topology )
    set l_orients [lindex [ analyze get_lipid_orients ] 1]

    set localorients [lindex [analyze lipid_orient_order all] 1]

    foreach or $localorients {
	set x [lindex $or 0]
	set y [lindex $or 1]

	# First figure out the absolute value of x and y components together
	set xyproj [expr sqrt($x*$x + $y*$y)]

	# Work out the bin to which this belongs
	if { $xyproj < $range } {
	    set thisbin [expr int(floor(($xyproj)/$binwidth))]
	    # Increment the relevant bin
	    lset av_localorients [namespace current]bin [expr [lindex $av_localorients $thisbin] + 1]
	}



    }

    incr av_localorients_i
    
}


