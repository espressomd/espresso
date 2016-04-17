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
# mbtools::utils -- 
#
# This package provides a set of routines for basic bookkeeping tasks
# that typically arise in a membrane simulation #
# Author: Ira
# 

namespace eval ::mbtools::utils {
    namespace export setup_outputdir
    namespace export read_startfile
    namespace export readcheckpoint
    namespace export read_topology
    namespace export set_topology
    namespace export set_bonded_interactions
    namespace export set_tab_interactions
    namespace export calculate_n_lipids
    namespace export init_random
    namespace export initialize_vmd
}

# ::mbtools::utils::setup_outputdir --
#
# This routine is designed to setup a directory for simulation
# output. It copies forcetables and the parameter file to the
# directory after creating it if necessary.  
#
# Arguments: outputdir The full pathname of the directory to be setup.
#                        At least the parent of this directory must
#                        exist
#
#            tabdir    The directory containing forcetables
#
#            tabnames  The list of all forcetables see
#                         set_tab_interactions
# 
#            paramsfile: The name of the parameter file 
#                        to be transferred to the directory
#
proc ::mbtools::utils::setup_outputdir { outputdir args } {

    # ---- Process Arguments ---------------------------# 
    # the ntabs option should not be necessary .. remove in a later
    # version
    set options {
	{paramsfile.arg  ""    "parameter file"}
	{tabdir.arg      ""    "forcetable directory"}
	{tabnames.arg    ""    "list of forcetable names"}
    }
    set usage "Usage: setup_outputdir \[paramsfile:tabdir:tabnames:startf:ntabs:] outputdir "
    # Strip off optional arguments and store in params
    array set params [::cmdline::getoptions args $options $usage]
    
    
    ::mmsg::send [namespace current]  "setting up $outputdir for $params(paramsfile)"
    # If <outputdir> doesn't exist then create it
    set errcode [catch { exec mkdir $outputdir  }]
    set ntabs [llength $params(tabnames)]
    
    # Copy forcetables to the current directory and to output directory
    for { set i 0 } { $i < $ntabs } { incr i } {
	set tablename [lindex $params(tabnames) $i ]
	set errcode [ catch { exec cp $params(tabdir)/$tablename [pwd]/ } ]	    
	if { $errcode } {
	    ::mmsg::warn [namespace current]  "couldn't transfer forcetable $params(tabdir)/$tablename to [pwd]"
	} else {
	    ::mmsg::send [namespace current]  "copied $params(tabdir)/$tablename to [pwd] "
	}

	set errcode [ catch { exec cp $params(tabdir)/$tablename $outputdir } ]
	if { $errcode } {
	    ::mmsg::warn [namespace current]  "couldn't transfer forcetable $params(tabdir)/$tablename to $outputdir"
	} else {
	    ::mmsg::send [namespace current]  "copied $params(tabdir)/$tablename to $outputdir "
	}
    }
    
    #Copy the paramsfile to the outputdir
    catch { exec cp $params(paramsfile) $outputdir }
    
    # Construct a directory for checkpoint backups inside outputdir
    catch { exec mkdir $outputdir/checkpoint_bak }    
}



# ::mbtools::utils::read_startfile --
#
# Read in configuration information from an existing config file
# 
#
proc ::mbtools::utils::read_startfile { file } {
    ::mmsg::send [namespace current]   "reading config: $file" nonewline 
    flush stdout
    set mark .gz
    set fname $file
    if { [ regexp $mark $fname ] } {
	set f [open "|gzip -cd $file" r]
    } else {
	set f [open "$file" r]
    }

    while { [ blockfile $f read auto ] != "eof" } {
	::mmsg::send [namespace current]  "." nonewline
	flush stdout
    }
    
    close $f
    ::mmsg::send [namespace current]  " done" 
}






# ::mbtools::utils::readcheckpoint --
#
#   A wrapper routine for reading checkpoints that checks for success
#
proc ::mbtools::utils::readcheckpoint { checkpointdir } {
    
    set result [ catch {  open "$checkpointdir/checkpoint.latest.gz" r  } ]
    if { $result } {
	::mmsg::warn [namespace current]  "no checkpoint named $checkpointdir/checkpoint.latest.gz"
	return 0
    }
    ::mmsg::send [namespace current] "reading Checkpoint $checkpointdir/checkpoint.latest.gz"
    #checkpoint_read "$checkpointdir/checkpoint" 0
    set in [open "|gzip -cd $checkpointdir/checkpoint.latest.gz" "r"]
    while { [blockfile $in read auto] != "eof" } {}
    close $in
    return 1
}

# ::mbtools::utils::read_topology --
#
#   Uses the blockfile read command to read topology information
#   contained in a file "$dir/file" into espresso.  The analyze
#   "topo_part_sync" command is then used to synchronise the molecule
#   information contained in topology and particle molecule id's
#
proc ::mbtools::utils::read_topology { file } {
    set f [open $file r]
    blockfile $f read topology
    analyze set "topo_part_sync"
}

# ::mbtools::utils::set_topology --
#
#  Takes a topology list "topo" and sets this into espresso. The
#   analyze "topo_part_sync" command is then used to synchronise the
#   molecule information contained in topology and particle molecule
#   id's
#
proc ::mbtools::utils::set_topology { topo } {
    eval analyze set $topo
    analyze set "topo_part_sync"
}


# ::mbtools::utils::set_bonded_interactions --
#
#  This routine uses the inter command of espresso to set all bonded
#  interactions.
#
proc ::mbtools::utils::set_bonded_interactions { bonded_parms } {
    foreach bondtype $bonded_parms {
	if { [catch {eval [concat inter $bondtype] } ] } {
	    mmsg::err [namespace current] "couldn't set interaction: [concat $int]"
	} else {	
	    mmsg::send [namespace current] "set interaction: $bondtype "
	}
    }
    return
}

# ::mbtools::utils::set_nb_interactions --
# 
# Set all the non-bonded interactions apart from tabulated ones, eg
# for constraints for instance.
#
proc ::mbtools::utils::set_nb_interactions { interactionlist } {
    foreach intertype $interactionlist {
	if { [catch { eval [concat inter  $intertype ] } ] } {
	    mmsg::err [namespace current] "could not set interaction: $intertype"
	}
	mmsg::send [namespace current] "set interaction: $intertype "
    }
    return
}

# ::puts_utilities::set_tab_interactions --
#
# DONT USE THIS ROUTINE ANYMORE  USE "set_nb_interactions" INSTEAD
#
# This routine sets up non_bonded interactions using tabulated
# potentials. 
#
# Arguments: tablenames A list containing the names of each tabulated
#             potential.  The list assumes an upper diagonal
#             structure.  For example if we have 3 particle types with
#             tabulated interactions then tablenames would be <t11 t12
#             t13 t22 t23 t33>
#           
# tabledir:  Main repository directory for force tables
#
# outputdir: Main directory for simulation output
#
# clean : Set this flag to yes if you want forcetables to be removed
#                from calling directory after they have been set
#
# Known bugs: Because table reading is performed by espresso in the
#                 current directory this can lead to problems when
#                 many scripts are started at once such as when using
#                 the opteron cluster.  In that case several instances
#                 of espresso may try to read the same file.  One
#                 solution to this would be to read the forcetables
#                 from the output directory which should be different
#                 for each instance of espresso.
#
proc ::mbtools::utils::set_tab_interactions { tablenames tabledir outputdir { clean "no" } } {

    variable maxnbinteractiontype

    # Find out how many forcetables we have	
    set ntabs [llength $tablenames]    
    

    for { set i 0 } { $i < $ntabs } {incr i} {
	# Copy the forcetables from tabledir to the current directory
	# where they will be read by espresso
	set errcode [catch { exec cp $tabledir/[lindex $tablenames $i] . } ]
	if { $errcode } {
	    ::mmsg::err [namespace current]  "couldn't copy table to [pwd]"
	}
	# Copy the forcetables from the current directory to the
	# output directory so that a record of the forcetables used is
	# kept with the output data
	set errcode [catch { exec cp [lindex $tablenames $i] $outputdir } ]
	if { $errcode } {
	    ::mmsg::err [namespace current]  "couldn't copy table to $outputdir "
	}
    }

    # Work out the number of atom types from the number of tables in <tablenames>
    set natomtypes [expr int(floor( ( sqrt(1+8*$ntabs)-1)/2.0))]

    # Set the interactions using the espresso command "inter"
    set tabcount 0
    for { set i 0 } { $i < $natomtypes } { incr i } {
	for { set j $i } {$j < $natomtypes } { incr j } {
	    inter $i $j tabulated [lindex $tablenames $tabcount]
	    mmsg::send [namespace current] "set tabulated interaction $i $j [lindex $tablenames $tabcount]"
	    incr tabcount
	    incr maxnbinteractiontype
	}
    }

    if { $clean == "yes" } {
	#clean up files from current directory
	for { set i 0 } { $i < [llength $tablenames] } {incr i} {
	    catch { exec rm [lindex $tablenames $i] }
	}
    }
    
}

    
# ::mbtools::utils::init_random --
#
# Initialize the random number generators on each processor 
#
proc ::mbtools::utils::init_random { n_procs } { 
    
    set c [expr [clock seconds] - 1068130800]    
    #    set c  1068130800    
    
    ::mmsg::send [namespace current]  "Setting the random seed to clock in seconds: $c "
    set seedbase $c
    for { set i 1 } { $i < $n_procs } { incr i } {
	lappend seedbase [expr $c + 2304*$i]
    }
    eval t_random seed $seedbase
    
    flush stdout
}
    
# ::mbtools::utils::initialize_vmd --
#
# Depending on the value of <flag> initialize vmd to one of two possible states:
#
#  interactive: VMD is started and a connection to espresso
#               established for immediate viewing of the current
#               espresso process. With some luck this might even work
#               sometimes!!!  If VMD doesn't get a proper connection
#               to espresso then it will crash.
#
#  offline: Just constructs the appropriate psf and
#               vmd_animation.script files and writes them to the
#               output directory so that pdb files generated with
#               writepdb can be viewed with vmd -e
#               vmd_animation.script
#
#  default: Any value other than those above for flag will just result
#               in vmd not being initialized.
#
proc ::mbtools::utils::initialize_vmd { flag outputdir ident  args } {
    # ---- Process Arguments ---------------------------# 
    set options {
	{extracommands.arg       ""     "additional stuff to be written to vmd_animation script"}
    }
    set usage "Usage: setup_outputdir \[paramsfile:tabdir:tabnames:startf:ntabs:] outputdir "
    # Strip off optional arguments and store in params
    array set params [::cmdline::getoptions args $options $usage]

    # set some defaults
    set filename "vmd"  
    set wait "0"
    set start "1" 
    
    ::mmsg::send [namespace current]  "initializing vmd to : " nonewline
    
    switch -regexp $flag  {
	"interactive" {
	    ::mmsg::send [namespace current]  "interactive" 
#	    writepsf_special "$filename.psf" 
	    writepst "$filename.psf"
	    writepdb "$filename.pdb"
	    for {set port 10000} { $port < 65000 } { incr port } {
		catch {imd connect $port} res
		if {$res == ""} break
	    }
	    set HOSTNAME [exec hostname]
	    set vmdout_file [open "vmd_start.script" "w"]
	    puts $vmdout_file "mol load psf $filename.psf pdb $filename.pdb"
	    puts $vmdout_file "rotate stop"
	    puts $vmdout_file "mol modstyle 0 0 lines"
	    # 1.800000 0.300000 8.000000 6.000000"
	    puts $vmdout_file "mol modcolor 0 0 SegName"
	    puts $vmdout_file "imd connect $HOSTNAME $port"
	    puts $vmdout_file "imd transfer 1"
	    puts $vmdout_file "imd keep 1"
	    close $vmdout_file
	    if { $start == 0 } {
		::mmsg::send [namespace current]  "Start VMD in the same directory on the machine you with :"
		::mmsg::send [namespace current]  "vmd -e vmd_start.script &"
		imd listen $wait
	    } else {
		exec vmd -e vmd_start.script &
	    }
	}
	"offline" {
	    ::mmsg::send [namespace current]  "offline"
	    variable firstconfignum
	    
#	    writepsf_special "$outputdir/$ident.vmd.psf" 
	    writepsf "$outputdir/$ident.vmd.psf"
	    
	    set vmd_file [open "$outputdir/vmd_animation.script" "w"]
	    puts $vmd_file "loadseries $ident.vmd 1 $firstconfignum"
	    puts $vmd_file "rotate stop"
	    puts $vmd_file "mol modstyle 0 0 lines"
	    puts $vmd_file "mol modcolor 0 0 SegName"
	    puts $vmd_file "logfile vmd.log"
	    puts $vmd_file "logfile off"
	    foreach command $params(extracommands) {
		puts $vmd_file $command
	    }
	    close $vmd_file
	    
	    # And for the warmup too
#	    writepsf_special "$outputdir/warm.vmd.psf" 
	    writepsf "$outputdir/warm.vmd.psf" 
	    set vmd_file [open "$outputdir/warm_animation.script" "w"]
	    puts $vmd_file "loadseries warm.vmd 1 0"
	    puts $vmd_file "rotate stop"
	    puts $vmd_file "mol modstyle 0 0 lines"
	    puts $vmd_file "mol modcolor 0 0 SegName"
	    puts $vmd_file "logfile vmd.log"
	    puts $vmd_file "logfile off"
	    foreach command $params(extracommands) {
		puts $vmd_file $command
	    }
	    close $vmd_file
	}


	"default" { 
	    ::mmsg::send [namespace current]  "default"
	    #Do nothing 
	}
    }
}






#-- Redundant routines -------------------------------#
# No longer used 

# ::mbtools::utils::writepsf_special --
#
# This will write a psf file that makes all the head beads the same
# regardless of whether they are really different.  Actually I would
# be carefull and probably never use this
#
#
proc ::mbtools::utils::writepsf_special { file {N_P -1} {MPC -1} {N_CI -1} {N_pS -1} {N_nS -1} } {
    set f [open $file "w"]
    puts $f "PSF"
    puts $f [format "%8d !NATOM" [setmd n_part]]
	# write atoms and create bondlist
    set cnt 1
    set mp [setmd max_part]

    
    set bondlist {}
    set bondcnt 0
    if { $N_P==-1 } {
	for {set p 0} { $p <= $mp } { incr p } {
	    set tp [part $p p t]
	    if { $p < $mp } {
		set tp1 [part [expr $p + 1] p t]
		set tp2 [part [expr $p + 2] p t]
	    }
	    if {  $tp == 0  &&  $tp1 == 2  } {
		set tp 3
	    }
	    if {  $tp == 0  &&  $tp2 == 2  } {
		set tp 3
	    }
	    
	    if { $tp != "na" } {
		set l [format "%8d T%03d %4d UNX  FE   FE"  $cnt $tp $p]
		if {[string length $l] < 50} {
		    set pad [string repeat " " [expr 50 - [string length $l]]]
		    set l "$l$pad"
		}
		puts $f $l
		set count($p) $cnt
		set bonds [part $p p b]
		foreach bb $bonds {
		    foreach b $bb {  # use 'set b [lindex $bb 0]' instead if you don't wanna have more than the 1st bond
			if {[llength $b ] == 2} { incr bondcnt; lappend bondlist "{$cnt [lindex $b 1]}" }
		    }
		}
	    }
	    incr cnt
	}
    } else {
	# include topology information in the data stream
	set ids 0
	set j1 [eval list        0         [expr $N_P*$MPC]       [expr $N_P*$MPC+$N_CI]       [expr $N_P*$MPC+$N_CI+$N_pS]       ]
	set j2 [eval list [expr $N_P*$MPC] [expr $N_P*$MPC+$N_CI] [expr $N_P*$MPC+$N_CI+$N_pS] [expr $N_P*$MPC+$N_CI+$N_pS+$N_nS] ]
	foreach ja $j1 je $j2 {
	    for {set p $ja} {$p < $je} {incr p} {
		set tp [part $p p t]
		if { $tp != "na" } {
		    set l [format "%8d T%03d %4d %03d  FE   FE"  $cnt $tp $p $ids]
		    if {[string length $l] < 50} {
			set pad [string repeat " " [expr 50 - [string length $l]]]
			set l "$l$pad"
		    }
		    puts $f $l
		    set count($p) $cnt
		    set bonds [part $p p b]
		    foreach bb $bonds {
			foreach b $bb {  # use 'set b [lindex $bb 0]' instead if you don't wanna have more than the 1st bond
			    if {[llength $b ] == 2} { incr bondcnt; lappend bondlist "{$cnt [lindex $b 1]}" }
			}
		    }
		}
		incr cnt; if { ($p < $N_P*$MPC-1) && !(($p+1) % $MPC) } { incr ids }
	    }
	    incr ids
	}
    }
    #  write bonds
    puts $f [format "%8d !NBOND                                      " $bondcnt]
    set bondlinecnt 0
    foreach b $bondlist {
	set b [lindex $b 0]
	if {[llength $b] == 2} {
	    incr bondcnt
	    eval "set p1 [lindex $b 0]"
	    eval "set p2 \$count([lindex $b 1])"
	    puts -nonewline $f [format "%8d%8d" $p1 $p2]
	    incr bondlinecnt
	    if {$bondlinecnt == 4} {
		puts $f ""
		set bondlinecnt 0
	    }
	}
    }
    close $f
}





# ::mbtools::utils::set_std_topology --
#    
# This routine is used when a starting config file is specified but no
# topology information is given.  In that case as a last resort we
# simply assume a homogenous bilayer ie standard topology.
#
proc ::mbtools::utils::set_std_topology { n_parts beads_per_lipid } {
    ::mmsg::send [namespace current] "assuming flat bilayer with one lipid type"
    set n_lipids_monolayer [expr $n_parts/($beads_per_lipid*2.0)]

    for {set i 0} { $i < $n_lipids_monolayer } {incr i} {
	set mol1 0
	set mol2 1
	for { set b 0 } { $b < $beads_per_lipid } {incr b } {
	    set partnum [expr $i*$beads_per_lipid*2 + $b ]
	    lappend mol1 $partnum
	}
	for { set b $beads_per_lipid } { $b < [expr $beads_per_lipid*2] } { incr b } {
	    set partnum [expr $i*$beads_per_lipid*2 + $b ]
	    lappend mol2 $partnum
	}
	lappend topo $mol1
	lappend topo $mol2
    }
    eval analyze set $topo
    analyze set "topo_part_sync"
}

    
# No longer used ?
proc ::mbtools::utils::prepare_vmd_series { outdir name  } {
    writepsf "$outdir/$name.vmd.psf" 
    
    set vmd_file [open "$outdir/vmd_animation.script" "w"]
    puts $vmd_file "loadseries $outdir/$name.vmd 1"
    puts $vmd_file "rotate stop"
    puts $vmd_file "mol modstyle 0 0 lines"
    puts $vmd_file "mol modcolor 0 0 SegName"
    puts $vmd_file "logfile $outdir/vmd.log"
    puts $vmd_file "logfile off"
    close $vmd_file
}

proc ::mbtools::utils::probe_nparts { f } {
    blockfile $f read variable
    blockfile $f read variable
    
    blockfile $f read interactions
    blockfile $f read particles  
    blockfile $f read bonds
    
    set n_lipids "[expr [ setmd n_part]/3]"
    return $n_lipids
}


proc ::mbtools::utils::cleanup_tabfiles { tablenames tabledir outputdir } {

    
    for { set i 0 } { $i < [llength $tablenames] } {incr i} {
	set errcode [catch { exec cp [lindex $tablenames $i] $outputdir } ]
	if { $errcode } {
	    ::mmsg::err [namespace current]  "Error copying table to $outputdir : $errcode"
	}
    }
    #clean up files
    for { set i 0 } { $i < [llength $tablenames] } {incr i} {
	catch { exec rm [lindex $tablenames $i] }
    }
}



