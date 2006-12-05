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
		    # default radius is 1.0
		    lset list $type 1.0
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
		if { $type <= [llength $name_list] } then {
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
	    return "atom $to $desc" 
	} else {
	    return "atom $from:$to $desc"
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

    set n_part [setmd n_part]

    # get the maximal type
    set type_max 0
    for { set pid 1 } { $pid < $n_part } { incr pid } {
	set type [part $pid print type]
	if { $type > $type_max } then { set type_max $type }
    }
    foreach {i j} $radius {
	if { $i > $type_max } then { set type_max $i }
    }
    foreach {i j} $typedesc {
	if { $i > $type_max } then { set type_max $i }
    }

    set typedesclist [get_typedesc_list $type_max $typedesc]
    set radiuslist [get_radius_list $type_max $radius]

    ##################################################
    # OUTPUT
    ##################################################
    # Print unitcell data
    if { $short } then { 
	puts $file "u [setmd box_l]"
    } else { 
	puts $file "unitcell [setmd box_l]"
    }

    # Print atom data

    set from 0
    set to 0
    set prev_type [part 0 print type]
    for { set pid 1 } { $pid < $n_part } { incr pid } {
	set type [part $pid print type]
	if { $prev_type == $type } then { 
	    incr to 
	} else {
	    # output from $from to $to
	    puts $file [get_atom_record $from $to $prev_type]

	    set from $pid
	    set to $pid
	    set prev_type $type
	}
    }
    puts $file [get_atom_record $from $to $prev_type]

    # Print bond data
    for { set from 0 } { $from < $n_part } { incr from } {
	set bonds [lindex [part $from print bond] 0]
	for { set i 0 } { $i < [llength $bonds] } { incr i } {
	    set to [lindex $bonds $i 1]

	    if { $short } then {
		puts $file "b $from:$to"
	    } else {
		puts $file "bond $from:$to"
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
	    default { 
		puts "unknown option to writevcf: $arg"
		return 
	    }
	}
    }

    # write the data
    set n_part [setmd n_part]

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
	    puts $file "u [setmd box_l]"
	} else { 
	    puts $file "unitcell [setmd box_l]"
	}
    }

    # output particle data
    if { $pids eq "all" } then {
	for { set pid 0 } { $pid < $n_part } { incr pid } {
	    if { $folded } then {
		puts $file [part $pid print folded_pos]
	    } else {
		puts $file [part $pid print pos]
	    }
	}
    } else {
	foreach pid $pids {
	    if { $folded } then {
		puts $file "$pid [part $pid print folded_pos]"
	    } else {
		puts $file "$pid [part $pid print pos]"
	    }
	}
    }

    if { ! $short } then { puts $file "" }
}
