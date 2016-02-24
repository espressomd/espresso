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
# ::mbtools::analysis::analyze_box_len --
#
# Extract the box dimensions from espresso
# 

namespace eval ::mbtools::analysis {}

namespace eval ::mbtools::analysis::boxl {
    variable av_boxl { 0 0 0 }
    variable av_boxl_i 0
    variable f_tvsbl
    variable verbose

    namespace export setup_boxl
    namespace export analyze_boxl
    namespace export printav_boxl
    namespace export resetav_boxl
}

proc ::mbtools::analysis::boxl::resetav_boxl { } {
    variable av_boxl 
    variable av_boxl_i 
    set av_boxl { 0.0 0.0 0.0 }
    set av_boxl_i 0
}

proc ::mbtools::analysis::boxl::printav_boxl { } {
    variable av_boxl_i
    variable av_boxl
    variable f_tvsbl
    global ::mbtools::analysis::time
    
    if { $av_boxl_i > 0 } {
	set avblx [expr [lindex $av_boxl 0]/($av_boxl_i*1.0)]
	set avbly [expr [lindex $av_boxl 1]/($av_boxl_i*1.0)]
	set avblz [expr [lindex $av_boxl 2]/($av_boxl_i*1.0)]
	puts $f_tvsbl "$time $avblx $avbly $avblz"
    } else {
	::mmsg::warn [namespace current] "can't print average box length"
	flush stdout
    }
    flush $f_tvsbl
}

proc ::mbtools::analysis::boxl::setup_boxl { args } {
    global  ::mbtools::analysis::outputdir
    variable f_tvsbl
    global ::mbtools::analysis::iotype
    global ::mbtools::analysis::suffix
    variable verbose
    ::mmsg::debug [namespace current] "opening $outputdir/time_vs_boxl$suffix "
    
     set options {
	{verbose "print out lots of stuff" }
    }
    set usage "Usage: setup_boxl verbose "
    array set params [::cmdline::getoptions args $options $usage]

    set verbose $params(verbose)

    if { [file exists "$outputdir/time_vs_boxl$suffix"] } {
	set newfile 0
    } else { 
	set newfile 1
    }
    set f_tvsbl [open "$outputdir/time_vs_boxl$suffix" $iotype]
    if { $newfile || $iotype == "w"} {
	puts $f_tvsbl "\# Time boxx boxy boxz"
    }
    

}

proc ::mbtools::analysis::boxl::analyze_boxl {  } {
    ::mmsg::send [namespace current] "analyzing box l"
    variable av_boxl_i
    variable av_boxl
    variable verbose

    set inst_boxl [setmd box_l]
    lset av_boxl 0 [expr [lindex $av_boxl 0] + [lindex $inst_boxl 0] ]
    lset av_boxl 1 [expr [lindex $av_boxl 1] + [lindex $inst_boxl 1] ]
    lset av_boxl 2 [expr [lindex $av_boxl 2] + [lindex $inst_boxl 2] ]

    incr av_boxl_i

    if { $verbose } {
	set avblx [expr [lindex $av_boxl 0]/($av_boxl_i*1.0)]
	set avbly [expr [lindex $av_boxl 1]/($av_boxl_i*1.0)]
	set avblz [expr [lindex $av_boxl 2]/($av_boxl_i*1.0)]
	::mmsg::send [namespace current]  "L: [lindex $inst_boxl 0] [lindex $inst_boxl 1] [lindex $inst_boxl 2] :: <L> $avblx $avbly $avblz"
	flush stdout
    }
    ::mmsg::debug [namespace current] "done"
}
