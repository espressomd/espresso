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
#-----------------------------------------------------------#
#
# Example script for using the data of an existing simulation.
# It runs the analysis routines out of the checkpoints that
# are stored in the simulation's folder.
#
# Implements settings contained in a parameter file which should be
# given as an argument.
#
# Note: You need to have consistent names between the configuration
# file name, the folder name which contains the checkpoints, and
# the checkpoint names. Also, it will obviously crash if you intend
# to change the initial conditions of the simulation! This does
# NOT recalculate the simulation, it merely runs any analysis routine
# you want. The only freedom you have is in changing the analysis
# flags, options, etc.
# 
# Enjoy!
#
# -----------------------------------------------------------#

# Get the name of the current namespace
set this [namespace current]


# --- ensure required packages are available  ---------------#
set result [package require ::mmsg]
 ::mmsg::send $this "loaded version [package require ::mbtools::analysis] of analysis"
 ::mmsg::send $this "loaded version [package require cmdline] of cmdline"
 ::mmsg::send $this "loaded version [package require ::mbtools::system_generation] of system_generation"

# Capture the child namespaces of analysis and system_generation so
# that we can explicitly allow messages from these.
set message_allowlist { :: ::mbtools::utils ::mbtools::system_generation ::mbtools::analysis }
set children [namespace children ::mbtools::analysis]
foreach child $children {
    lappend message_allowlist $child
}
set children [namespace children ::mbtools::system_generation]
foreach child $children {
    lappend message_allowlist $child
}

# Set the namespaces from which messages will be printed
catch { ::mmsg::setnamespaces $message_allowlist }




# ---- Process Command Line Args --------------------------- #

# Note that this can cause problems reading checkpoints.  The usage
# below seems to work ok.
set options {
    {n.arg      1    set the number of processors }
}
set usage "Usage: main.tcl n: paramsfile"
array set params [::cmdline::getoptions argv $options $usage]
set paramsfile [lindex $argv 0]

# Enable debugging messages
#::mmsg::enable debug

#----------- System Parameters ----------------------------#
# Set a bunch of default parameters.
set thermo Langevin
set warmup_temp 0
set warmsteps 0
set warmtimes 0
set free_warmsteps 0 
set free_warmtimes 0 
set startj 0    
set startk 0
set startmdtime 0
set npt off
set mgrid 8
set stray_cut_off 1000.0
set use_vmd "offline"
set tablenames ""
set trappedmols ""
::mmsg::send $this "using paramsfile: $paramsfile"
source $paramsfile




# ----------- Initialization ------------------ -----------#

# Attempt to read a checkpoint file - check there is at least one checkpoint!
set checkpointexists [ ::mbtools::utils::readcheckpoint $outputdir ]

# Set the starting time for this run ie override the value in checkpoint file
set starttime [clock seconds]

if { !$checkpointexists } {
    mmsg::err $this "Couldn't read any checkpoint."
} else {
    # Make sure the analysis_flags are taken from $paramsfile
    unset analysis_flags
    source $paramsfile
    # A checkpoint exists so all we need to do is reset the topology and setup analysis again
    ::mbtools::utils::read_topology "$outputdir/$topofile"
    ::mbtools::analysis::setup_analysis $analysis_flags -outputdir  $outputdir -g $mgrid -str $stray_cut_off

    # Yikes I hope this is right.  We want to make sure that we start
    # exactly from where the checkpoint was written
    set startj $j
    set startk [expr $k + 1]

}

#Main Integration                                          #
#----------------------------------------------------------#

set j $startj

set timingstart [clock clicks -milliseconds]

set filecounter 0
set fileswitch 0
while { !$fileswitch } {
    
    set fileswitch [ catch { open "$outputdir/$ident.[format %04d $filecounter].gz" r }]
    if { !$fileswitch } {
	mmsg::send $this "reading checkpoint $filecounter"

	# Easy way to get all the checkpoints available to be read.
	set pipe [open "$outputdir/checkpoint_tmp.chk" w ]
	puts $pipe "$outputdir/$ident.[format %04d $filecounter].gz"
	close $pipe

	checkpoint_read "$outputdir/checkpoint_tmp" 0

	# Again, source $paramsfile to have the newest variables
	source $paramsfile
	
	# Call all of the analyze routines that we specified when setting
	# up our analysis
	::mbtools::analysis::do_analysis
	::mbtools::analysis::print_averages
	
    }
    set filecounter [expr $filecounter + 1]

}

# terminate program
mmsg::send $this "\n\nfinished"
exit 1



