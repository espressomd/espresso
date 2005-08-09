#!/bin/sh
# \
    exec $ESPRESSO_SOURCE/Espresso $0 1 $*
# series.tcl [output]
#    files <file 0> <file 1> .....
#    pattern <prefix> <start> <stop> [<skip> [<pad>]]
#####################################################

set filename vmd
set argv [lrange $argv 1 end]

if {[lindex $argv 0]=="output"} {
    set argv [lrange $argv 1 end]
    set output 1
} else {set output 0}

switch [lindex $argv 0] {
    files {
	set argv [lrange $argv 1 end]
	if {$output==1} {
	    set i 0
	    foreach e $argv {
		part deleteall
		polyBlockRead "$e"
		if {$i == 0} {writepsf "$filename.psf"}
		writepdb "$filename[format %04d $i].pdb"
		incr i 
	    }
	    set vmdout_file [open "vmd_start.script" "w"]
	    puts $vmdout_file "logfile vmd.log"
	    puts $vmdout_file "rotate stop"
	    puts $vmdout_file "logfile off"
	    puts $vmdout_file "mol modstyle 0 0 CPK 1.800000 0.300000 8.000000 6.000000"
	    puts $vmdout_file "mol modcolor 0 0 SegName"
	    puts $vmdout_file "source $env(ESPRESSO_SCRIPTS)/vmd_plg.tcl"
	    puts $vmdout_file "loadseries $filename"
	    close $vmdout_file

	    exec vmd -e vmd_start.script &
	} else {
	    set i 0
	    foreach e $argv {
		part deleteall
		polyBlockRead "$e"
		if {$i == 0} {prepare_vmd_connection $filename 3000}
		imd positions
		incr i 
	    }
	   
	}
    }
    pattern {
	set prefix [lindex $argv 1]
	set start [lindex $argv 2]
	set stop [lindex $argv 3]
	if {[lindex $argv 4]!=0} {set skip [lindex $argv 4]} else {set skip 1}
	if {[lindex $argv 5]!=0} {set pad [lindex $argv 5]} else {set pad 4}
	if {$output==1} {
	    set j 0
	    for {set i $start} {$i <= $stop} {incr i} {
		part deleteall
		polyBlockRead "$prefix[format %0${pad}d $i]"
		if {$i == $start} {writepsf "$filename.psf"}
		writepdb "$filename[format %04d $j].pdb"
		incr j
	    }
	    set vmdout_file [open "vmd_start.script" "w"]
	    puts $vmdout_file "logfile vmd.log"
	    puts $vmdout_file "rotate stop"
	    puts $vmdout_file "logfile off"
	    puts $vmdout_file "mol modstyle 0 0 CPK 1.800000 0.300000 8.000000 6.000000"
	    puts $vmdout_file "mol modcolor 0 0 SegName"
	    puts $vmdout_file "source $env(ESPRESSO_SCRIPTS)/vmd_plg.tcl"
	    puts $vmdout_file "loadseries $filename"
	    close $vmdout_file

	    exec vmd -e vmd_start.script &
	} else {
	    set j 0
	    for {set i $start} {$i <= $stop} {incr i} {
		part deleteall
		polyBlockRead "$prefix[format %0${pad}d $i]"
		if {$i == $start} {prepare_vmd_connection $filename 3000}
		imd positions
		incr j
	    }
	}
    }
    default {puts "ERROR: switch must be either files or pattern;\nAborting...";exit}
}


