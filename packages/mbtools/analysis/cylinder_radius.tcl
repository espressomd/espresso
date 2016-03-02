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
# ::mbtools::analysis::analyze_cylinder_radius --
#
# Calculate the radius of a cylinder
#   

namespace eval ::mbtools::analysis {}

namespace eval ::mbtools::analysis::cylinder_radius {
    variable R_middle 0
    variable R_inner 0
    variable R_outer 0

    variable verbose

    namespace export printav_cylinder_radius
    namespace export setup_cylinder_radius
    namespace export analyze_cylinder_radius
    namespace export resetav_cylinder_radius
}

proc ::mbtools::analysis::cylinder_radius::resetav_cylinder_radius { } {
    variable R_middle 0
    variable R_inner 0
    variable R_outer 0
}

proc ::mbtools::analysis::cylinder_radius::printav_cylinder_radius { } {
    variable R_middle
    variable R_inner
    variable R_outer
    variable f_tvsen
    global ::mbtools::analysis::time
    puts -nonewline $f_tvsen "$time $R_middle $R_inner $R_outer"
    puts $f_tvsen ""
    flush $f_tvsen
}

proc ::mbtools::analysis::cylinder_radius::setup_cylinder_radius { args } {
    global ::mbtools::analysis::outputdir
    global ::mbtools::analysis::suffix
    global ::mbtools::analysis::iotype
    variable f_tvsen
    variable verbose

    set options {
	{verbose "print out lots of stuff" }
    }
    set usage "Usage: setup_cylinder_radius verbose "
    array set params [::cmdline::getoptions args $options $usage]

    set verbose $params(verbose)

    if { [file exists "$outputdir/time_vs_cylinder_radius$suffix"] } {
	set newfile 0
    } else { 
	set newfile 1
    }
    ::mmsg::debug  [namespace current]  "opening $outputdir/time_vs_cylinder_radius$suffix "
    set f_tvsen [open "$outputdir/time_vs_cylinder_radius$suffix" $iotype ]

    if { $newfile || $iotype == "w"} {
	puts $f_tvsen "\# Components of the cylinder_radius "
	puts $f_tvsen "\# Time R_middle R_inner R_outer"
    }
    flush $f_tvsen
}


proc ::mbtools::analysis::cylinder_radius::analyze_cylinder_radius {  } {
    ::mmsg::send [namespace current] "analyzing cylinder_radius"
    variable R_middle
    variable R_inner
    variable R_outer
    variable verbose

    set box [setmd box_l]

    set topology [analyze set]

    set com {0. 0.}
    set num_mol 0
    # Assume cylinder is built along the z-direction
    # Calculate center of mass of the cylinder in the x-y plane.    
    foreach mol $topology {
	set headp [lindex $mol 1]
	set tb [lindex $mol end]
	
	# NONfolded particle position is needed
	set hpos [part $headp print pos]

	# Check the neighborhood of this particle and see if it is a stray particle
	set nb [analyze nbhood $headp 3]
	if { [llength $nb] > 5 } {
	    # Not stray particle
	    lset com 0 [expr [lindex $com 0]+[lindex $hpos 0]]
	    lset com 1 [expr [lindex $com 1]+[lindex $hpos 1]]
	    incr num_mol
	}
    }
    lset com 0 [expr [lindex $com 0]/(1.*$num_mol)]
    lset com 1 [expr [lindex $com 1]/(1.*$num_mol)]

    if { $verbose } {
	::mmsg::send [namespace current] "COM: [lindex $com 0] [lindex $com 1]"
    }
    set num_outer 0
    # Now calculate two radii: inner and outer by looking at the orientation of a lipid
    foreach mol $topology {
	set headp [lindex $mol 1]
	set tb [lindex $mol end]

	set hpos [part $headp print pos]
	# Only look at x,y coordinates
	set hpos [lrange $hpos 0 end-1]
	set tpos [part $tb    print pos]
	# Only look at x,y coordinates
	set tpos [lrange $tpos 0 end-1]

	# Define vector CoM to head bead
	set com_h_vec [::mbtools::utils::min_vec $hpos $com]

	# Define vector tail bead to head bead

	set t_h_vec [::mbtools::utils::min_vec $hpos $tpos]

	# dot product between the two vectors
	set com_h_dot_t_h [vecdot_product $com_h_vec $t_h_vec]
	
	if {$com_h_dot_t_h > 0} {
	    incr num_outer
	    # same orientation -- it's the outer leaflet
	    set R_outer [expr $R_outer + [veclen $com_h_vec]]
	} else {
	    # antiparallel orientation -- it's the inner leaflet
	    set R_inner [expr $R_inner + [veclen $com_h_vec]]
	}
    }
	
    set R_outer [expr $R_outer/(1.*$num_outer)]
    set R_inner [expr $R_inner/(1.*($num_mol-$num_outer))]
    set R_middle [expr ($R_inner+$R_outer)/2.]

    if { $verbose } {
	::mmsg::send [namespace current] "cylinder_radius: R_middle: $R_middle -- R_inner: $R_inner -- R_outer: $R_outer"
    }
    
    ::mmsg::debug [namespace current] "done"

}


