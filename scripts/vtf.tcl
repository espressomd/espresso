#
# Copyright (C) 2012,2013 The ESPResSo project
# Copyright (C) 2006,2007,2008,2009,2010,2011 Olaf Lenz
#  
# This file is part of ESPResSo.
#  
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 
#  
#############################################################
#                                                           #
# vtf.tcl                                                   #
# =======                                                   #
#                                                           #
# Functions that allow writing VTF files.                   #
#                                                           #
#############################################################

# Write the structure information of the current system to the given file
# OPTIONS:
#   short|verbose [verbose]
#   radius [{}]
#   typedesc [{}]
proc writevsf { file args } {
    proc flatten_indexed_list { n default_value indexed_list } {
	set flatlist {}
	# create a list of length n
	for { set i 0 } { $i < $n } { incr i } {
	    set flatlist [lappend flatlist $default_value]
	}

	foreach {index value} $indexed_list {
	    lset flatlist $index $value
	}

	return $flatlist
    }

    proc get_radius_list { type_max i_list } {
	set list [flatten_indexed_list [expr $type_max + 1] "auto" $i_list]
	for { set type 0 } { $type < [llength $list] } { incr type } {
	    if { [lindex $list $type] eq "auto" } then {
		set l [split [inter $type $type]]
		if {[lindex $l 2] eq "lennard-jones" } then {
		    # guess the radius from the self-lj-interaction
		    set lj_sigma [lindex $l 4]
		    set lj_offset [lindex $l 7]
		    set r [expr 0.5*($lj_sigma+$lj_offset)]
		    lset list $type $r
		} else {
		    # default radius is 0.5
		    lset list $type 0.5
		}
	    }
	}
	return $list
    }

    proc get_typedesc_list { type_max i_list } {
	set name_list { "O" "N" "S" "H" "C" "Z" }
	set list [flatten_indexed_list [expr $type_max + 1] "default" $i_list]
	for { set type 0 } { $type < [llength $list] } { incr type } {
	    if { [lindex $list $type] eq "default" } then {
		if { $type < [llength $name_list] } then {
		    lset list $type "name [lindex $name_list $type] type $type"
		} else {
		    lset list $type "name X type $type"
		}
	    }
	}
	return $list
    }

    # returns an atom record for the atoms $from to $to of type $type
    proc get_atom_record { from to type } {
	upvar 1 radiuslist radiuslist typedesclist typedesclist
	set desc "radius [lindex $radiuslist $type] [lindex $typedesclist $type]"
	if { $to == $from } then { 
	    return "atom [vtfpid $to] $desc" 
	} else {
	    return "atom [vtfpid $from]:[vtfpid $to] $desc"
	}
    }

    ##################################################
    # PROCEDURE CODE
    ##################################################
    set typedesc {}
    set radius {}
    set short 0

    # Parse options
    for { set argnum 0 } { $argnum < [llength $args] } { incr argnum } {
	set arg [ lindex $args $argnum ]
	set val [ lindex $args [expr $argnum+1]]
	switch -- $arg {
	    "typedesc" { 
		set typedesc $val
		incr argnum 
	    }
	    "radius" { 
		set radius $val
		incr argnum 
	    }
	    "verbose" { set short 0 }
	    "short" { set short 1 }
	    default { 
		puts "unknown option to writevsf: $arg"
		return 
	    }
	}
    }

    # cleanup particle mapping
    global vtf_pid
    array unset vtf_pid

    # Print unitcell data
    if { $short } then { 
	puts $file "p [setmd box_l]"
    } else { 
	puts $file "pbc [setmd box_l]"
    }

    set max_pid [setmd max_part]

    # compute the maximal type
    set max_type [setmd n_part_types]
    foreach {i j} $radius {
	if { $i > $max_type } then { set max_type $i }
    }
    foreach {i j} $typedesc {
	if { $i > $max_type } then { set max_type $i }
    }

    set typedesclist [get_typedesc_list $max_type $typedesc]
    set radiuslist [get_radius_list $max_type $radius]

    # collect combinations of charge and type
    set uses_electrostatics [has_feature "ELECTROSTATICS"]
    set qs ""
    for { set pid 0 } { $pid <= $max_pid } { incr pid } {
        set type [part $pid print type]
        if { $uses_electrostatics } then { 
            set q [part $pid print q]
            set qs "q$q"
        }
        set this "t$type$qs"
        if { ! [info exists combination($this)] } then {
            set combination($this) $pid
            set desc "radius [lindex $radiuslist $type]"
            set desc "$desc [lindex $typedesclist $type]"
            if { $uses_electrostatics } then { set desc "$desc q $q" }
            set combination("desc-$this") $desc
        } else {
            lappend combination($this) $pid
        }
    }

    # loop over the combinations and create atom record
    foreach c [array names combination "t*"] {
        set pids $combination($c)
        set desc $combination("desc-$c")
        set to [vtfpid [lindex $pids 0]]
        set aids ""
        foreach pid $pids {
            set vpid [vtfpid $pid]
            if { [expr $vpid-1] != $to } then {
                if { [info exists from] } then {
                    # print out records from $from to $to
                    if { $from == $to } then {
                        lappend aids "$to"
                    } else {
                        lappend aids "$from:$to"
                    }
                }
                set from $vpid
            }
            set to $vpid
        }
        if { $from == $to } then {
            lappend aids "$to"
        } else {
            lappend aids "$from:$to"
        }
        unset from
        # let's group atom ranges, so that there are no more than 8 per line
        # This fixes a problem with the vmd plugin, and makes the vtf more 
        # readable anyway.
        set start 0
        set maxlen 8
        set ll  [llength [lrange $aids $start end ]]
        while { $ll >= $maxlen } {
              puts $file "atom [join [lrange $aids $start [expr $start + $maxlen -1]] ,] $desc"
              incr start $maxlen
              set ll  [llength [lrange $aids $start end ]]
        }
        if { $start < [llength $aids ] } { 
          puts $file "atom [join [lrange $aids  $start  end] ,] $desc"
        }
    }
    
    # Print bond data
    for { set from 0 } { $from <= $max_pid } { incr from } {
	if { [part $from] != "na" } then {
	    set bonds [lindex [part $from print bond] 0]
	    for { set i 0 } { $i < [llength $bonds] } { incr i } {
		set to [lindex $bonds $i 1]
		
		if { $short } then {
		    puts $file "b [vtfpid $from]:[vtfpid $to]"
		} else {
		    puts $file "bond [vtfpid $from]:[vtfpid $to]"
		}
	    }
	}
    }

    if { ! $short } then { puts $file "" }
}

# Write the coordinates of the system to the file
# OPTIONS:
#   short|verbose - [verbose]
#   folded|absolute - write the coordinates folded or not [absolute]
#   pids <pids> - write only the coordinates of these particles [all]
proc writevcf { file args } {
    # set defaults
    set folded 0
    set short 0
    set pids "all"
    set userdata 0

    # Parse options
    for { set argnum 0 } { $argnum < [llength $args] } { incr argnum } {
	set arg [ lindex $args $argnum ]
	set val [ lindex $args [expr $argnum + 1]]
	switch -- $arg {
	    "folded" { set folded 1 }
	    "absolute" { set folded 0 }
	    "verbose" { set short 0 }
	    "short" { set short 1 }
	    "pids" { set pids $val; incr argnum }
	    "userdata" { set userdata $val; incr argnum }
	    default { 
		puts "unknown option to writevcf: $arg"
		return 
	    }
	}
    }

    # write the data
    set max_pid [setmd max_part]

    if { $pids eq "all" } then {
	if { $short } \
	    then { puts $file "t" } \
	    else { puts $file "timestep ordered" }
    } else {
	if { $short } \
	    then { puts $file "i" } \
	    else { puts $file "timestep indexed" }
    }

    # unitcell
    if { [lindex [integrate] 0 1] eq "npt_isotropic" } then {
	if { $short } then { 
	    puts $file "p [setmd box_l]"
	} else { 
	    puts $file "pbc [setmd box_l]"
	}
    }

    # output particle data
    if { $pids eq "all" } then {
	for { set pid 0 } { $pid <= $max_pid } { incr pid } {
	    if {[part $pid] != "na"} then {
		if { $folded } then {
		    puts -nonewline $file [part $pid print folded_pos]
		} else {
		    puts -nonewline $file [part $pid print pos]
		}
		if { [llength $userdata] } then {
		    puts $file [ lindex $userdata $pid ]
		} else {
		    puts $file ""
		}
	    }
	}
    } else {
	foreach pid $pids {
	    if {[part $pid] != "na"} then {
		if { $folded } then {
		    puts -nonewline $file "[vtfpid $pid] [part $pid print folded_pos]"
		} else {
		    puts -nonewline $file "[vtfpid $pid] [part $pid print pos]"
		}
		if { [llength $userdata] } then {
		    puts $file [ lindex $userdata $pid ]
		} else {
		    puts $file ""
		}
	    }
	}
    }

    if { ! $short } then { puts $file "" }
}

# get the VMD pid of a given ESPResSo-PID
proc vtfpid { pid } {
    global vtf_pid
    if {! [array exists vtf_pid]} then {
	# compute the mapping from Espresso PIDS to VMD PIDS
	set i 0
	set max_pid [setmd max_part]
	for { set epid 0 } { $epid <= $max_pid } { incr epid } {
	    if {[part $epid] == "na"} then {
		set vtf_pid($epid) "na"
	    } else {
		set vtf_pid($epid) $i
		incr i
	    }
	}
    }
    return $vtf_pid($pid)
}

