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
# ::mbtools::analysis::analyze_pressure --
#
#  Calculate the total pressure of the system and break it into components
#

namespace eval ::mbtools::analysis {}

namespace eval ::mbtools::analysis::pressure {
    variable av_pressure { 0 0 0 0 0 0 }
    variable av_pressure_i 0
    variable f_tvsp
    variable verbose 0

    namespace export setup_pressure
    namespace export analyze_pressure
    namespace export printav_pressure
    namespace export resetav_pressure

}

proc ::mbtools::analysis::pressure::resetav_pressure { } {
    variable av_pressure 
    variable av_pressure_i 
    set av_pressure {0.0 0.0 0.0 0.0 0.0 0.0 }
    set av_pressure_i 0
}

proc ::mbtools::analysis::pressure::printav_pressure { } {
    variable av_pressure_i
    variable av_pressure
    variable f_tvsp
    global ::mbtools::analysis::time

    if { $av_pressure_i > 0 } {
	puts -nonewline $f_tvsp "$time "
	for { set v 0 } { $v < [llength $av_pressure] } {incr v} {
	    puts -nonewline $f_tvsp "[expr [lindex $av_pressure $v]/($av_pressure_i*1.0)] "
	    
	}
	puts $f_tvsp " "
    } else {
	::mmsg::warn [namespace current] "can't print average pressure"
	flush stdout		       
    }
    flush $f_tvsp
}


proc ::mbtools::analysis::pressure::setup_pressure { args } {
    global ::mbtools::analysis::outputdir
    global ::mbtools::analysis::suffix
    global ::mbtools::analysis::iotype
    variable f_tvsp
    variable verbose

    ::mmsg::debug [namespace current] "Opening $outputdir/time_vs_pressure$suffix "
    
    set options {
	{verbose "print out lots of stuff" }
    }
    set usage "Usage: setup_pressure verbose "
    array set params [::cmdline::getoptions args $options $usage]

    set verbose $params(verbose)

    if { [file exists "$outputdir/time_vs_pressure$suffix"] } {
	set newfile 0
    } else { 
	set newfile 1
    }
    
    set f_tvsp [open "$outputdir/time_vs_pressure$suffix" $iotype]
    if { $newfile || $iotype == "w"} {
	puts $f_tvsp "\# Components of the total pressure"
	puts $f_tvsp "\# Time p_inst(NPT only) total ideal fene harmonic nonbonded "
    }
    
}

proc ::mbtools::analysis::pressure::analyze_pressure {  } {
    ::mmsg::send [namespace current] "analyzing pressure"
    variable av_pressure
    variable av_pressure_i
    variable verbose
    set p_inst [setmd npt_p_inst_av]
    set p_tmp [analyze pressure]
    lset av_pressure 0 [expr [lindex $av_pressure 0] + $p_inst]
    for { set i 0 } { $i < [llength $p_tmp ] } { incr i } {
	set tmp [lindex $p_tmp $i]
	set ntmp [llength $tmp]	
	if { [ regexp "pressure" $tmp ] } {
	    lset av_pressure 1 [expr [lindex $av_pressure 1] + [lindex $tmp 1]]
	}
	if { [ regexp "ideal" $tmp ] } {
	    lset av_pressure 2 [expr [lindex $av_pressure 2] + [lindex $tmp 1]]
	}
	if { [ regexp "FENE" $tmp ] } {
	    lset av_pressure 3 [expr [lindex $av_pressure 3] + [lindex $tmp 2]]
	}
	if { [ regexp "HARMONIC" $tmp ] } {
	    lset av_pressure 4 [expr [lindex $av_pressure 4] + [lindex $tmp 2]]
	}
	if { [ regexp "nonbonded" $tmp ] } {
	    lset av_pressure 5 [expr [lindex $av_pressure 5] + [lindex $tmp 3]]
	}
	
    }
    incr av_pressure_i
    if { $verbose } {
	::mmsg::send [namespace current] "p_inst [lindex $av_pressure 0] : pressure [lindex $av_pressure 1] : ideal [lindex $av_pressure 2] : FENE [lindex $av_pressure 3] : HARMONIC [lindex $av_pressure 4] : nonbonded [lindex $av_pressure 5]"	    
    }
    ::mmsg::debug [namespace current] "done"
}
