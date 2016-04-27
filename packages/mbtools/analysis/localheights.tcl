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
# Procedures for determining the distribution of local relative lipid heights
#

namespace eval ::mbtools::analysis {}

namespace eval ::mbtools::analysis::localheights {
    variable av_localheights
    variable av_localheights_i 0
    variable nbins
    variable range
    variable rcatch
    variable verbose
    namespace export setup_localheights
    namespace export analyze_localheights
    namespace export printav_localheights
    namespace export resetav_localheights
}

proc foldpart { pos box } {
    for { set i 0 } { $i < 3 } { incr i } {
	set x  [lindex $pos $i]
	set x [expr $x - floor($x/[lindex $box $i])*[lindex $box $i] ]
	lset pos $i $x
    }
    return $pos
}

#
# Returns a list of the molecule ids of all lipids within a certain radius
# 
# l_orients: A list of all orient flags for all molecules obtained
#           from a call to analyze get_lipid_orients
# 
# orient : Flag for top or bottom leaflet (see analyze
#        get_lipid_orients for definitions of this flag
#
#
# headp : Atom ID of central bead
# 
# rcatch: Catch radius in xy plane
#
# molid : The molecule id of the central molecule

proc nbhoodlipids { l_orients orient headp rcatch molid } {
    
    # Search for the neighbours around this head bead
    set nblist [analyze nbhood planar 1 1 0 $headp $rcatch  ]

    # Turn the bead neighbourlist into a molecule neighbour list
    set molnblist -1 
    foreach bead $nblist {
	set thismol [part $bead print mol]
	# If the molecule id is not yet in our mol list then add it
	# We also check if the lipid is in the same monolayer and only include it if it is
	
	#		puts "$orient [lindex $l_orients $thismol]"
	if { [lsearch $molnblist $thismol] == -1 && $thismol != $molid && [lindex $l_orients $thismol] == $orient } {
	    lappend molnblist $thismol
	    if { [lindex $molnblist 0 ] == -1 } {
		set molnblist $thismol
	    }		
	}
    }	    
    return $molnblist
}

proc ::mbtools::analysis::localheights::resetav_localheights  { } {
    # Do nothing

}

proc ::mbtools::analysis::localheights::printav_localheights {  } {
    variable av_localheights
    variable av_localheights_i
    variable f_localheights
    variable nbins
    variable range
    global ::mbtools::analysis::topology
    global ::mbtools::analysis::outputdir
    global ::mbtools::analysis::suffix

    if { $av_localheights_i > 0 } {
	
	set f_localheights [open "$outputdir/av_localh$suffix" "w" ]
	# Write a header to the file
	puts $f_localheights "\# local height distribution of lipids "
	puts $f_localheights ""

	
	set lhbinwidth [expr $range/($nbins*1.0)]
	
	# Write out the local height distribution to a file
	for { set bin 0 } { $bin < [llength $av_localheights] } { incr bin } {
	    set currbin [lindex $av_localheights $bin]
	    puts $f_localheights "[expr $bin*$lhbinwidth - $range/2.0] [expr $currbin/(1.0*[llength $topology])]"
	}		    
	close $f_localheights
    }
}

proc ::mbtools::analysis::localheights::setup_localheights { args } {
    variable nbins
    variable range
    variable verbose
    variable rcatch

    set options {
	{verbose "print out lots of stuff" }
	{range.arg "1.0" "Range for localorients calculation" }
	{nbins.arg "100" "Number of bins used for analysis" }
	{rcatch.arg "1.9" "Catch radius for finding the neighbourhood of lipids" }
    }
    set usage "Usage: setup_localheights verbose:nbins:range:rcatch "
    array set params [::cmdline::getoptions args $options $usage]
    
    set verbose $params(verbose)
    set range $params(range)
    set nbins $params(nbins)
    set rcatch $params(rcatch)

    set av_localheights 0
    unset av_localheights
    #Initialize av_localheights
    for { set bn 0 } { $bn < $nbins } { incr bn } {
	lappend av_localheights 0
    }
    
}

proc ::mbtools::analysis::localheights::calc_localheights {} {
    global ::mbtools::analysis::topology
    variable rcatch
    variable verbose
    set cmx 0
    set cmy 0
    set cmz 0

    set box [setmd box_l]

    # Internal copy of l_orients ( same molecule order as topology )
    set l_orients [lindex [ analyze get_lipid_orients ] 1]

    # Check that hdiffs is uninitialized
    set hdiffs 0
    unset hdiffs

    set molid 0
    foreach mol $topology {

	set headp [lindex $mol 1]
	set tb [lindex $mol end ]
	# For some old simulations where 0 was a special hb type we need this
	if { [part $headp print type] != 0 && [part $tb print type] == 0 } {
	    set headp $tb
	}
	set orient [lindex $l_orients $molid]

	# Do not perform any analysis for stray beads
	if { $orient < 2 } {

#	    puts "[part $headp print pos]"
	    set hpos [foldpart [part $headp print pos] $box ]

#	    puts "$molid $headp"


	    set tries 0
	    set maxtries 100
	    set molnblist ""
	    while { [llength $molnblist] != 6 && $tries < $maxtries } {
		if { [llength $molnblist] < 6 } {
		    set rcatch [expr $rcatch + 0.01 ]
		} else {
		    set rcatch [expr $rcatch - 0.01 ]
		}

		if { $rcatch <= 0 } {
		    ::mmsg::err [namespace current] "catch radius went to zero when trying to find lipid neighbours"
		}

		set molnblist [nbhoodlipids $l_orients $orient $headp $rcatch $molid ]
#		puts "[llength $molnblist] $tries $rcatch"
		incr tries

	    }

	    # Use the molnblist to work out height differences
	    foreach nbmol $molnblist {
		if { $nbmol != -1 } {
		    set hb [lindex [lindex $topology $nbmol] 1 ]
		    set tb [lindex [lindex $topology $nbmol] end ]
#		    puts "$nbmol $hb $tb"
		    if { [ catch { set dum $hb } ] || [ catch { set dum $tb } ] } {
			::mmsg::err [namespace current] "no hb"
		    }
		
		    # For some old simulations where 0 was a special hb type we need this
		    if { [part $hb print type] != 0 && [part $tb print type] == 0 } {
			set hb $tb
		    }
		    set hbpos [foldpart [part $hb print pos] $box]
		    lappend hdiffs [expr [lindex $hpos 2] - [lindex $hbpos 2]]
		    #		puts "[lindex $topology $nbmol] [part $hb print type] $hpos $hbpos $hb"
		}
	    }
	    
#	    puts "$molnblist $hdiffs"
	    
	}
	incr molid
    }

    return $hdiffs

}

# ::mbtools::analysis::analyze_localheights --
# 
# Calculates the relative heights of all nearest neighbours of each lipid and bins these into a histogram
#
proc ::mbtools::analysis::localheights::analyze_localheights {  } {
    variable av_localheights
    variable av_localheights_i
    variable nbins
    variable range

    set nbins $nbins
    set binwidth [expr $range/($nbins*1.0)]

    ::mmsg::send [namespace current] "analyzing local heights"

    set localheights [calc_localheights]

    # Bin up the localheights
    foreach height $localheights {

	set thisbin [expr int(floor(($height+$range/(2.0))/$binwidth))]
#	puts "$height $thisbin"
	if { $thisbin > 0 && $thisbin < $nbins } {
#	    puts "$thisbin [llength $av_localheights]"
	    lset av_localheights $thisbin [expr [lindex $av_localheights $thisbin] + 1 ]
	}
    }

    incr av_localheights_i
    
}


