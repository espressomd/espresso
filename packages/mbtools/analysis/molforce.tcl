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
# ::mbtools::analysis::analyze_molforce --
#
# Finds the forces applied by a trap upon a molecule
# Espresso returns the total force applied by the trap summed over all time steps since
# the  "analyze mol force $mol" command was last called along with the number of time steps
# To find the average force the total force is divided by the number of timesteps

namespace eval ::mbtools::analysis {}

namespace eval ::mbtools::analysis::molforce {
    variable av_molforce
    variable f_tvsmf
    variable verbose

    namespace export setup_molforce
    namespace export analyze_molforce
    namespace export printav_molforce
    namespace export resetav_molforce
}

proc ::mbtools::analysis::molforce::resetav_molforce { } {
    variable av_molforce
    variable mollist

    set av_molforce ""
    foreach mol $mollist {
       lappend av_molforce { 0.0 0.0 0.0 0 }
    }
}

proc ::mbtools::analysis::molforce::printav_molforce { } {
    variable av_molforce
    variable f_tvsmf
    global ::mbtools::analysis::time
    variable mollist
    
    set i 0;
    puts $f_tvsmf "$time " nonewline
    foreach mol $mollist {
	set this_av_molforce [lindex $av_molforce $i]
	if {[lindex $this_av_molforce 3] > 0} {
	    # before printing the force we divided it by the number of time steps
	    set avmfx [expr [lindex $this_av_molforce 0]/([lindex $this_av_molforce 3]*1.0)]
	    set avmfy [expr [lindex $this_av_molforce 1]/([lindex $this_av_molforce 3]*1.0)]
	    set avmfz [expr [lindex $this_av_molforce 2]/([lindex $this_av_molforce 3]*1.0)]
	    puts $f_tvsmf "$avmfx $avmfy $avmfz " nonewline
	} else {
	    ::mmsg::warn [namespace current] "can't print average molforce"
	    flush stdout
	}
	incr i
    } 
    puts $f_tvsmf ""
    flush $f_tvsmf
}

proc ::mbtools::analysis::molforce::setup_molforce { args } {
    global  ::mbtools::analysis::outputdir
    variable f_tvsmf
    global ::mbtools::analysis::iotype
    global ::mbtools::analysis::suffix
    variable verbose
    variable av_molforce
    variable mollist

    ::mmsg::debug [namespace current] "opening $outputdir/time_vs_molforce$suffix "

    set options {
	{verbose "print out lots of stuff" }
	{mollist.arg "" "list of molecules analyzed"}
    }
    set usage "Usage: setup_molforce verbose "
    array set params [::cmdline::getoptions args $options $usage]

    set verbose $params(verbose)
    set mollist $params(mollist)

    if { [file exists "$outputdir/time_vs_molforce$suffix"] } {
	set newfile 0
    } else { 
	set newfile 1
    }
    set f_tvsmf [open "$outputdir/time_vs_molforce$suffix" $iotype]
    if { $newfile || $iotype == "w"} {
	puts $f_tvsmf "\#      " nonewline
	foreach mol $mollist {
	    puts $f_tvsmf "$mol           " nonewline
	}
	puts $f_tvsmf ""
	puts $f_tvsmf "\# Time " nonewline
	foreach mol $mollist {
	    puts $f_tvsmf "x   y   z   " nonewline
	}
	puts $f_tvsmf ""
    }
    foreach mol $mollist {
	lappend av_molforce { 0.0 0.0 0.0 0 }
    }
}

proc ::mbtools::analysis::molforce::analyze_molforce {  } {
    ::mmsg::send [namespace current] "analyzing molforce"
    variable av_molforce
    variable verbose
    variable mollist
    set av_molforce_new ""
    set i 0
    foreach mol $mollist {
	set inst_molforce [analyze mol force $mol]
	set this_av_molforce [lindex $av_molforce $i]
# This first three entries in $inst_molforce are the force components.  The fourth entry is the number of time step over which the force was summed.
	set this_av_molforce "[expr [lindex $this_av_molforce 0] + [lindex $inst_molforce 0] ]  \
                              [expr [lindex $this_av_molforce 1] + [lindex $inst_molforce 1] ]  \
                              [expr [lindex $this_av_molforce 2] + [lindex $inst_molforce 2] ]  \
                              [expr [lindex $this_av_molforce 3] + [lindex $inst_molforce 3] ]"
	lappend av_molforce_new $this_av_molforce
	if { $verbose } {
	    set ifx [expr [lindex $inst_molforce 0]/([lindex $inst_molforce 3]*1.0)]
	    set ify [expr [lindex $inst_molforce 1]/([lindex $inst_molforce 3]*1.0)]
	    set ifz [expr [lindex $inst_molforce 2]/([lindex $inst_molforce 3]*1.0)]
	    set avmfx [expr [lindex $this_av_molforce 0]/([lindex $this_av_molforce 3]*1.0)]
	    set avmfy [expr [lindex $this_av_molforce 1]/([lindex $this_av_molforce 3]*1.0)]
	    set avmfz [expr [lindex $this_av_molforce 2]/([lindex $this_av_molforce 3]*1.0)]
	    ::mmsg::send [namespace current]  "For mol $mol - L: $ifx $ify $ifz from [lindex $inst_molforce 3] steps :: <L> $avmfx $avmfy $avmfz from [lindex $this_av_molforce 3] steps"
	    flush stdout
	}
	incr i
    }   
    set av_molforce $av_molforce_new

    ::mmsg::debug [namespace current] "done"
}
