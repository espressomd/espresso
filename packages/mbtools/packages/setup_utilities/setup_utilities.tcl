# setup_utilities -- 
#
# This package provides a set of routines for basic bookkeeping tasks
# that typically arise in a membrane simulation #
# Author: Ira
# 

#package require mmsg 0.1.0
package provide setup_utilities 0.1.0

namespace eval ::setup_utilities {
    variable verbosity 0
    variable firstconfignum 0

    # The identity number for written out warmup configurations
    variable warmcfg 0

    variable maxnbinteractiontype 0

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

source [file join [file dirname [info script]] warmup.tcl]





# ::setup_utilities::setup_outputdir --
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
#            startf    An optional configuration file from which to
#                        extract forcetable info using
#                        "get_forcetablenames"
#
proc ::setup_utilities::setup_outputdir { outputdir args } {

    # ---- Process Arguments ---------------------------# 
    # the ntabs option should not be necessary .. remove in a later
    # version
    set options {
	{paramsfile.arg  ""    "parameter file"}
	{tabdir.arg      ""    "forcetable directory"}
	{tabnames.arg    ""    "list of forcetable names"}
	{startf.arg      ""    "an initial configuration file"}
	{ntabs.arg       6     "number of forcetables"}
    }
    set usage "Usage: setup_outputdir \[paramsfile:tabdir:tabnames:startf:ntabs:] outputdir "
    # Strip off optional arguments and store in params
    array set params [::cmdline::getoptions args $options $usage]
    
    
    ::mmsg::send [namespace current]  "setting up $outputdir for $params(paramsfile)"
    # If <outputdir> doesn't exist then create it
    set errcode [catch { exec mkdir $outputdir  }]
    
    
    if { $params(startf) != "" } {
	mmsg::send [namespace current] "you specified a starting config $params(startf)"
	# If a starting file is specified then copy it to the output directory 
	set errcode [ catch { exec cp $params(startf) $outputdir } ]
	# Extract the tables from the starting file : see
	# get_forcetablenames for possible problems with this
	set tablenames [ get_forcetablenames $fpath $ntabs ] 
    } else {
	# If no starting file is specified then we should ensure that
	# tablenames contains enough elements
	set ntabs [llength $params(tabnames)]
	if { $ntabs == 0 } { 
	    ::mmsg::err [namespace current]  "<tablenames> not set" 
	}
    }
    
    # Copy forcetables to the current directory
    for { set i 0 } { $i < $ntabs } { incr i } {
	set tablename [lindex $params(tabnames) $i ]
	set errcode [ catch { exec cp $params(tabdir)/$tablename [pwd]/ } ]	    
	if { $errcode } {
	    ::mmsg::warn [namespace current]  "couldn't transfer forcetable $params(tabdir)/$tablename to [pwd]"
	} else {
	    ::mmsg::send [namespace current]  "copied $params(tabdir)/$tablename to [pwd] "
	}
    }
    
    #Copy the paramsfile to the outputdir
    catch { exec cp $params(paramsfile) $outputdir }
    
    # Construct a directory for checkpoint backups inside outputdir
    catch { exec mkdir $outputdir/checkpoint_bak }    
}



# ::setup_utilities::read_startfile --
#
# Read in configuration information from an existing config file
# 
#
proc ::setup_utilities::read_startfile { file } {
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


# ::setup_utilities::get_forcetablenames --
#
# Given an existing output file <filepath> retrieve the list of forcetables from it.
#
# Todo: The <ntabs> argument should be removed since it can be
#         obtained from the file itself anyway
#
proc ::setup_utilities::get_forcetablenames { filepath ntabs } {
    set mark .gz
    if { [ regexp $mark $filepath ] } {
	set f [open "|gzip -cd $filepath" r]
    } else {
	set f [open "$filepath" r]
    }
    
    set data [read $f] 
    set num [lsearch -regexp $data .tab ] 
    set data [lindex $data $num]
    set num [lsearch -regexp $data .tab ] 
    
    for { set i 0 } { $i < $ntabs } { incr i } {
	set tmpdata [lindex $data [expr $num + $i]]
	lappend tables [lindex $tmpdata 3]
    }
    return $tables
}



# ::setup_utilities::readcheckpoint --
#
#   A wrapper routine for reading checkpoints that checks for success
#
proc ::setup_utilities::readcheckpoint { checkpointdir } {
    
    set result [ catch {  open "$checkpointdir/checkpoint.latest.gz" r  } ]
    if { $result } {
	::mmsg::warn [namespace current]  "no checkpoint named $checkpointdir/checkpoint.latest.gz"
	return 0
    }
    ::mmsg::send [namespace current] "reading Checkpoint $checkpointdir/checkpoint.latest.gz"
    checkpoint_read "$checkpointdir/checkpoint" 0
    return 1
}

# ::setup_utilities::read_topology --
#
#   Uses the blockfile read command to read topology information
#   contained in a file "$dir/file" into espresso.  The analyze
#   "topo_part_sync" command is then used to synchronise the molecule
#   information contained in topology and particle molecule id's
#
proc ::setup_utilities::read_topology { file } {
    set f [open $file r]
    blockfile $f read topology
    analyze set "topo_part_sync"
}

# ::setup_utilities::set_topology --
#
#  Takes a topology list "topo" and sets this into espresso. The
#   analyze "topo_part_sync" command is then used to synchronise the
#   molecule information contained in topology and particle molecule
#   id's
#
proc ::setup_utilities::set_topology { topo } {
    eval analyze set $topo
    analyze set "topo_part_sync"
}


# ::setup_utilities::set_topology --
#
#  This routine uses the inter command of espresso to set all bonded
#  interactions.
#
# Arguments: 
#
#  bonded_params A tcl list of the form: Interactionname param1 param2
#                   ... where the parameters should correspond exactly
#                   to that required by espresso.  Possible values for
#                   Interaction name are;
#
#   FENE <k_fene> <r_fene>
#   Harm_bend <k_harmonic> <r_harmonic>
#   Harm <k_harmonic <r_harmonic>
#
#  Note that Harm_bend and Harm are actually just different names for
#  the same interaction type.  Harm_bend should be applied between
#  atoms separated by a bridging atom that is itself FENE bonded to
#  the others.  This results in an effective bending stiffness
#  interaction that is pairwise.  
#
#
# Note that at present this just sets up a completely standard set of
# bonded interactions which is used for all the systems in the
# simulation. Ultimately we would like to be able to make this
# different for every lipid.  To do that we would need to modify the
# way lipids are placed.
#
proc ::setup_utilities::set_bonded_interactions { bonded_parms } {
    set num 0
    foreach val $bonded_parms {
	switch [lindex $val 0] {
	    "FENE" {
		inter $num   FENE  [lindex $val 1]  [lindex $val 2]
	    }
	    "Harm_bend" {
		inter $num harmonic [lindex $val 1] [lindex $val 2]
	    }
	    "Harm" {
		inter $num harmonic [lindex $val 1] [lindex $val 2]
	    }
	}
	incr num
    }
}

# ::setup_utilities::set_nb_interactions --
# 
# Set all the non-bonded interactions apart from tabulated ones, eg
# for constraints for instance.
#
proc ::setup_utilities::set_nb_interactions { interactionlist } {
    for { set i 0 } {$i <  [llength $interactionlist]} { incr i } {
	set int [lindex $interactionlist $i]
	if { [catch {eval [concat inter  $int]} ] } {
	    mmsg::err [namespace current] "couldn't set interaction: [concat $int]"
	} else {	
	    mmsg::send [namespace current] "set interaction: [inter [lindex $int 0] [lindex $int 1]]"
	}
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
proc ::setup_utilities::set_tab_interactions { tablenames tabledir outputdir { clean "no" } } {

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




# ::setup_utilities::calculate_n_lipids --
#
# This routine calculates the number of lipids we would require to
# create a square sheet of length <len> given a certain area fraction
# and beadsize.
#   
proc ::setup_utilities::calculate_n_lipids { fraction  beadsize len avflag } {
    set lx [lindex $len 0]
    set ly [lindex $len 1]
    set lz [lindex $len 2]
    
    set area [expr $lx*$ly]
    set volume [expr $lx*$ly*$lz]
    
    if { $avflag == "-af" } {
	set n_lipids_monolayer [expr floor($area*$fraction/(3.1415*0.25*$beadsize*$beadsize))]
	set n_lipids [expr $n_lipids_monolayer*2]	
    } elseif { $avflag == "-vf" } {
	::mmsg::err [namespace current]  "vf Not yet implemented "
    } else { 
	::mmsg::err [namespace current]  "calculate_n_lipids No such flag: $avflag" 
    }
    
    return $n_lipids
}

    
# ::setup_utilities::init_random --
#
# Initialize the random number generators on each processor up to a maximum of <n_procs> = 8
#
proc ::setup_utilities::init_random { n_procs } { 
    
    set c [expr [clock seconds] - 1068130800]    
    #    set c  1068130800    
    
    ::mmsg::send [namespace current]  "Setting the random seed to clock in seconds: $c "
    
    
    if { $n_procs == 1} { t_random seed $c 
    } elseif { $n_procs == 2 } { t_random seed $c [expr $c + 100] 
    } elseif { $n_procs == 3 } { t_random seed $c [expr $c + 100] [expr $c + 200]  } elseif { $n_procs == 4 } { t_random seed $c [expr $c + 100] [expr $c + 200] [expr $c + 300] 
    }     elseif { $n_procs == 5 } { t_random seed $c [expr $c + 100] [expr $c + 200] [expr $c + 300] [expr $c + 400]
    }     elseif { $n_procs == 6 } { t_random seed $c [expr $c + 100] [expr $c + 200] [expr $c + 300] [expr $c + 400] [expr $c + 500]
    }    elseif { $n_procs == 7 } { t_random seed $c [expr $c + 100] [expr $c + 200] [expr $c + 300] [expr $c + 400] [expr $c + 500] [expr $c + 600]
    }    elseif { $n_procs == 8 } { t_random seed $c [expr $c + 100] [expr $c + 200] [expr $c + 300] [expr $c + 400] [expr $c + 500] [expr $c + 600] [expr $c + 700]
    }    else { 
	::mmsg::err [namespace current]  " Max n_procs exceeded " 
    }
    flush stdout
}
    
# ::setup_utilities::initialize_vmd --
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
proc ::setup_utilities::initialize_vmd { flag outputdir ident  args } {
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

# ::setup_utilities::writepsf_special --
#
# This will write a psf file that makes all the head beads the same
# regardless of whether they are really different.  Actually I would
# be carefull and probably never use this
#
#
proc ::setup_utilities::writepsf_special { file {N_P -1} {MPC -1} {N_CI -1} {N_pS -1} {N_nS -1} } {
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


# ::setup_utilities::set_std_topology --
#    
# This routine is used when a starting config file is specified but no
# topology information is given.  In that case as a last resort we
# simply assume a homogenous bilayer ie standard topology.
#
proc ::setup_utilities::set_std_topology { n_parts beads_per_lipid } {
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


#-- Redundant routines -------------------------------#
    
# No longer used ?
proc ::setup_utilities::prepare_vmd_series { outdir name  } {
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

proc ::setup_utilities::probe_nparts { f } {
    blockfile $f read variable
    blockfile $f read variable
    
    blockfile $f read interactions
    blockfile $f read particles  
    blockfile $f read bonds
    
    set n_lipids "[expr [ setmd n_part]/3]"
    return $n_lipids
}


proc ::setup_utilities::cleanup_tabfiles { tablenames tabledir outputdir } {

    
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