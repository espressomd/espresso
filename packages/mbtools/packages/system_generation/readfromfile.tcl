# 
#  Routines for generating bilayer systems
#


namespace eval system_generation {}




#::system_generation::readfromfile
#
# Read particle data from a file
#
#
# This proc does a crude read of the file and its topology and sets
# the values into espresso.  Note that there are two topologies here,
# the topology that is read in from the file and the topology that was
# generated automatically by the make_simple_topology command.  Some
# simple checks are made to ensure that these two correspond to each
# other.
#
# Bonds and bond types from the file will be set.  The routine returns
# a topology which is a corrected topology read from the topofile.
# Generally speaking the data that is read in should correspond in
# terms of molecule types, bonds, bond_types, particle number and
# topology to the values which would be generated from your input
# parameters.
#
proc ::system_generation::readfromfile {topo boxl filename topofile args } {
    ::mmsg::send [namespace current] "reading particle data from $filename "

    set options {
	{ignore.arg     ""   "particle properties to be ignored"  }
    }
    set usage "Usage: create_bilayer topo boxl \[bondl:uniform:]"
    array set params [::cmdline::getoptions args $options $usage]
    

    # Numerical tolerance for differences in box length between file
    # and simulation values
    set bxtol 0.000001

    # Check to see if the file is gzipped and unzip it before opening
    # if it is
    set mark .gz
    if { [ regexp $mark $filename ] } {
	set fd [open "|gzip -cd $filename" r]
    } else {
	set fd [open "$filename" r]
    }
    set ft [open "$topofile" ]

    # Read the raw data
    set raw [read $fd]
    set fraw [lindex [read $ft] 0]
    set ftopo [lrange $fraw 1 end]
    # Extract the particle and bonding data
    foreach item $raw {	
	if { [lsearch $item "particles" ] != -1 } {
	    set pdatakey [lindex $item 1]
	    set pdata [lrange $item 2 end]
	}
	if { [lsearch $item "bonds" ] != -1 } {
	    set bdata [lrange $item 1 end]
	}	
	if { [lsearch -regexp $item "box_l" ] != -1 } {
	    set fboxl [lrange [lindex $item 1] 1 end] 
	    set boxmismatch 0
	    for { set d 0 } { $d < 3 } { incr d } {
		set bxdif [expr [lindex $fboxl $d] - [lindex $boxl $d]]
		if { [expr $bxdif*$bxdif] > $bxtol } {
		    set boxmismatch 1
		}
	    }
	    if { $boxmismatch } {
		mmsg::warn [namespace current] "box dimensions $fboxl in file do not match the values $boxl set in this simulation"
	    }
	}
    }

    # Now check to see if the number of particles in the file matches
    # that in the topology 
    set natoms 0
    foreach mol $topo {
	set natoms [expr $natoms + [expr [llength $mol ] - 1] ]
    }
    if { $natoms != [llength $pdata] } {
	mmsg::err [namespace current] "[llength $pdata] atoms were specified in the file $filename but $natoms were specified in topology"
    }

    # We allow for the case that the topo part id's start from a value
    # other than zero and this will min value is added to the partids
    set minidtopo [minpartid $topo]
#    puts "[minmoltype $topo]  [minmoltype $ftopo]"

    set typeoffset [expr [minmoltype $topo] - [minmoltype $ftopo]]
    if {$typeoffset < 0 } { 
	set typeoffset [expr -$typeoffset] 
    }
    


    # Now place the atoms using the order designated in the file
    # itself and adjust the read-in topology to fit that of topo

    # Note that particle info is assumed to be 
    set m 0
    foreach mol $ftopo {
	for { set b 1 } { $b < [llength $mol] } { incr b } {
	    set bid [lindex $mol $b]
	    set pt [lindex $pdata $bid]
	    

	    # Construct the part command to add this particle
	    set k 0
	    set vl 0


	    while { $k < [llength $pdatakey] }  {
		switch -exact [lindex $pdatakey $k] {
		    "id" {
			if { $bid != [lindex $pt $vl] } {
			    mmsg::err [namespace current] "non consecutive particle numbers in file"
			}
			set val [expr [lindex $pt $vl] + $minidtopo]
			lappend pargs $val
			# Correct the pid in ftopo
			lset ftopo $m $b $val
			incr k
			incr vl
		    }
		    "pos" {
			# Not ignorable
			lappend pargs pos
			lappend pargs [lindex $pt $vl]
			incr vl
			lappend pargs [lindex $pt $vl]
			incr vl
			lappend pargs [lindex $pt $vl]
			incr vl
			incr k
		    }
		    "type" {
			# Not ignorable
			lappend pargs type
			set val [expr [lindex $pt $vl] + $typeoffset]
			lappend pargs $val
			incr vl
			incr k
		    }
		    "v" {
			if { [lsearch -exact $params(ignore) v] != -1 } {
			    incr vl 3
			    incr k
			} else {
			    lappend pargs v
			    lappend pargs [lindex $pt $vl]
			    incr vl
			    lappend pargs [lindex $pt $vl]
			    incr vl
			    lappend pargs [lindex $pt $vl]
			    incr vl
			    incr k
			}
		    }
		    "f" {
			if { [lsearch -exact $params(ignore) f] != -1 } {
			    incr vl 3
			    incr k
			} else {
			    lappend pargs f
			    lappend pargs [lindex $pt $vl]
			    incr vl
			    lappend pargs [lindex $pt $vl]
			    incr vl
			    lappend pargs [lindex $pt $vl]
			    incr vl
			    incr k
			}

		    }
		    "molecule" {
			# This gets set later so let's leave it
			incr k
			incr vl
		    }
		    "default" {

			# If unknown particle attributes exist in the
			# file just skip them (hopefully they are only
			# a single element)
			mmsg::debug [namespace current] "unrecognised particle attribute [lindex $pdatakey $k] assumed to have just 1 element"
			incr k
			incr vl
		    }
		}
	    }
	    eval [concat part $pargs]
	    unset pargs
	}
	incr m
    }
    

    # Now setup all the bonds from the file
    for {set bd 0 } { $bd < [llength $bdata] } { incr bd } {
	set thispart [lindex $bdata $bd 0]
	set thisbonds [lindex $bdata $bd 1] 
	set pid [expr $thispart + $minidtopo] 
	lappend bargs $pid
	lappend bargs bond
	for { set i 0 } { $i < [llength $thisbonds] } { incr i } {	    
	    set bondedto [lindex $thisbonds $i]
	    lset bondedto 1 [expr [lindex $bondedto 1] + $minidtopo]
	    eval [concat part $bargs $bondedto]
	}
	unset bargs
    }

    return $ftopo

}
