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
# Proceedures for determining the vertical density profile of a bilayer
#

namespace eval ::mbtools::analysis {}

namespace eval ::mbtools::analysis::radial_density_profile {
    variable av_densities
    variable av_tdensities
    variable av_densities_i 0
    variable verbose
    variable xbins
    variable ybins
    variable tbins
    variable yrange
    variable xrange
    variable beadtypes

    namespace export printav_radial_density_profile
    namespace export setup_radial_density_profile
    namespace export analyze_radial_density_profile
    namespace export resetav_radial_density_profile
}

proc ::mbtools::analysis::radial_density_profile::resetav_radial_density_profile { } {
    # Do nothing because we want to average this over the entire simulation
    
}

proc ::mbtools::analysis::radial_density_profile::printav_radial_density_profile { } {
    variable f_densityprof
    variable av_densities_i
    variable beadtypes
    variable av_densities
    variable av_tdensities
    variable xbins
    variable ybins
    variable tbins
    variable yrange
    variable xrange
    global ::mbtools::analysis::outputdir
    global ::mbtools::analysis::suffix

    if {  $av_densities_i > 0 } {	
	for { set bt 0 } { $bt < [llength $beadtypes] } { incr bt } {
	    set f_densityprof [open "$outputdir/av_rmap$bt$suffix" "w" ]
	    for { set by 0 } { $by <  $ybins } { incr by } {
		for { set bx 0 } { $bx <  $xbins } { incr bx } {
		    puts -nonewline $f_densityprof " [lindex $av_densities $bt $bx $by ]"
		}
		puts $f_densityprof ""
	    }
	    close $f_densityprof
	    set f_densityprof [open "$outputdir/av_tprof$bt$suffix" "w" ]
	    for { set ti 0 } { $ti < $tbins } { incr ti } {
		puts $f_densityprof "[lindex $av_tdensities $bt $ti ]"
	    }
	    close $f_densityprof
	}
    }
}

proc ::mbtools::analysis::radial_density_profile::setup_radial_density_profile { args } {

    global ::mbtools::analysis::outputdir
    global ::mbtools::analysis::suffix
    global ::mbtools::analysis::iotype
    variable xbins
    variable ybins
    variable tbins
    variable yrange
    variable xrange
    variable beadtypes
    variable av_densities
    variable av_tdensities
    variable av_densities_i

    set options {
	{xbins.arg "50" "Number of bins used for horizontal axis" }
	{ybins.arg "50" "Number of bins used for vertical axis" }
	{tbins.arg "100" "Number of bins used for theta" }
	{xrange.arg "50.0" "horizontal range over which to calculate density profile"}
	{yrange.arg "50.0" "horizontal range over which to calculate density profile"}
	{beadtypes.arg "0" "Identity of beads to use for density profile analysis" }
    }
    set usage "Usage: setup_radial_density_profile xbins:ybins:xrange:yrange:beadtypes "
    array set params [::cmdline::getoptions args $options $usage]

    set beadtypes $params(beadtypes)
    set xbins $params(xbins)
    set ybins $params(ybins)
    set tbins $params(tbins)
    set xrange $params(xrange)
    set yrange $params(yrange)

    #Initialize av_densities
    set thisxbinlist 0.0
    set thisybinlist 0.0
    set thistbinlist 0.0
    unset thisxbinlist
    unset thisybinlist
    unset thistbinlist
    for { set bt 0 } { $bt < [llength $beadtypes] } { incr bt } {
	for { set bx 0 } { $bx <  $xbins } { incr bx } {
	    for { set by 0 } { $by <  $ybins } { incr by } {
		lappend thisybinlist 0.0
	    }
	    lappend thisxbinlist $thisybinlist
	    unset thisybinlist
	}
	lappend av_densities $thisxbinlist
	unset thisxbinlist
	
	for { set ti 0 } { $ti < $tbins } { incr ti } {
	    lappend thistbinlist 0.0
	}
	lappend av_tdensities $thistbinlist
	unset thistbinlist
    }
    

    mmsg::send [namespace current] "setup radial density profile for beads $beadtypes with bindimensions $xbins x $ybins and range $xrange x $yrange"

}

# ::mbtools::analysis::analyze_radial_density_profile --
# 
# Calculates the density of each lipid bead type as a function of
# vertical distance relative to the bilayer midplane
#
proc ::mbtools::analysis::radial_density_profile::analyze_radial_density_profile { } {
    variable av_densities
    variable av_tdensities
    variable av_densities_i
    variable xbins
    variable ybins
    variable tbins
    variable yrange
    variable xrange
    variable beadtypes


    ::mmsg::send [namespace current] "analyzing radial density profile " nonewline

    set center [analyze centermass 0 ]
    # ultimately this should be obtained from the principle axes of the system
    set axes [analyze find_principal_axis [lindex $beadtypes 0]  ]
    set axis [lindex $axes 1 1]
    ::mmsg::send [namespace current] "with axis $axis" 

    set alldensities [analyze radial_density_map $xbins $ybins  $xrange $yrange [lindex $axis 0] [lindex $axis 1] [lindex $axis 2] [lindex $center 0] [lindex $center 1] [lindex $center 2] $beadtypes $tbins ]

    # Select the actual data and skip the label
    set densities [lindex $alldensities 0 1]
    
    # Bin up the data for the density map
    for { set bt 0 } { $bt < [llength $beadtypes] } { incr bt } { 
	for { set bx 0 } { $bx < $xbins } { incr bx } {
	    for { set by 0 } { $by < $ybins } { incr by } {
		lset av_densities $bt $bx $by [ expr [lindex $densities $bt $bx $by] + [lindex $av_densities $bt $bx $by] ]
	    }
	}
    }
    set densities [lindex $alldensities 1 1]

    # Bin up the data for the theta profile 
    for { set bt 0 } { $bt < [llength $beadtypes] } { incr bt } { 
	for { set ti 0 } { $ti < $tbins } { incr ti } {
	    lset av_tdensities $bt $ti [expr [lindex $densities $bt $ti] + [lindex $av_tdensities $bt $ti]]
	}
    }

    incr av_densities_i

    
}


