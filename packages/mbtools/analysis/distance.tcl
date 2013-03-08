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
# ::mbtools::analysis::analyze_distance --
#
# Calculate the distance between two particles.  At present this
# calculates distance between two "middlebeads" on molecules.  In
# principle this routine should be generalized to calculate the
# distance between any two particles
# 
# Author Gregoria

namespace eval ::mbtools::analysis {}

namespace eval ::mbtools::analysis::distance {  
    variable f_tvsdist
    variable av_dist  0
    variable av_dist_i 0
    namespace export printav_distance
    namespace export setup_distance
    namespace export analyze_distance
    namespace export resetav_distance
}
proc ::mbtools::analysis::distance::resetav_distance {} {
    variable av_dist  
    variable av_dist_i 
    set av_dist 0.0
    set av_dist_i 0   
}

proc ::mbtools::analysis::distance::printav_distance { } {
    variable f_tvsdist
    variable av_dist  0
    variable av_dist_i 0
    global ::mbtools::analysis::time

    if {  $av_dist_i > 0 } {
	set avdistx [expr [lindex $av_dist 0]/($av_dist_i*1.0)]
	
	puts $f_tvsdist "$time $avdistx"
	#                   set av_boxl_i 0
    } else {
	::mmsg::warn [namespace current] "can't print average distance"
	flush stdout
    }
    flush $f_tvsdist
}

proc ::mbtools::analysis::distance::setup_distance { args } {
    global ::mbtools::analysis::outputdir
    global ::mbtools::analysis::suffix
    global ::mbtools::analysis::iotype

    ::mmsg::debug [namespace current] "opening $outputdir/time_vs_distance$suffix "
    
    set options {
	{verbose "print out lots of stuff" }
    }
    set usage "Usage: setup_distance verbose "
    array set params [::cmdline::getoptions args $options $usage]

    set verbose $params(verbose)
	
    if { [file exists "$outputdir/time_vs_distance$suffix"] } {
	set newfile 0
    } else {
	set newfile 1
    }
    set f_tvsdist [open "$outputdir/time_vs_distance$suffix" $iotype]
    if { $newfile || $iotype == "w"} {
	puts $f_tvsdist "\# Time Distance"
    }
}




proc ::mbtools::analysis::analyze_distance { } {

    ::mmsg::send [namespace current] "analyzing the distance between 2 particles"
    variable av_dist_i
    variable av_dist

    variable middlebead
    
    set middlebead [::system_generation::get_middlebead ]

    set inst_dist [bond_length [lindex $middlebead 0] [lindex $middlebead 1]] 

    lset av_dist 0 [expr [lindex $av_dist 0] + [lindex $inst_dist 0] ]

    #lset av_dist 1 [expr [lindex $av_dist 1] + [lindex $inst_dist 1] ]
    #lset av_dist 2 [expr [lindex $av_dist 2] + [lindex $inst_dist 2] ]

    incr av_dist_i

  
    if { $verbose } {
	set avdistx [expr [lindex $av_dist 0]/($av_dist_i*1.0)]

	#set avdisty [expr [lindex $av_dist 1]/($av_dist_i*1.0)]
	#set avdistz [expr [lindex $av_dist 2]/($av_dist_i*1.0)]



	::mmsg::send [namespace current]  "L: [lindex $inst_dist 0] :: <L> $avdistx"
	flush stdout
    }
    ::mmsg::debug [namespace current] "done"
}
