######################################
# particle support
######################################

proc blockfile_write_particles {channel write particles {info "id pos v q"} {range "0-end"}} {
    blockfile $channel write start "particles"
    if {![regexp "pos" $info]} {set info "pos $info"}
    if {![regexp "id" $info]} {set info "id $info"}
    if {[regexp "bonds" $info]} { error "bonding information cannot be written" }

    puts $channel "{$info} "
    foreach list $range {
	set list [split $list "-"]
	set len [llength $list]
	if {$len == 1} {
	    if {$list == "all"} {
		set start 0
		set end [setmd max_part]
	    } {
		if {$list == "end"} {set list [setmd max_part]}
		set start $list
		set end   $list
	    }
	} {
	    if {$len != 2} {
		error "list must have form \"<a>-<b> <c>-end...\""
		exit
	    }
	    set start [lindex $list 0]
	    if {$start == "end"} {set start [setmd max_part]}
	    set end [lindex $list 1]
	    if {$end == "end"} {set end [setmd max_part]}
	}
	if {![string is integer $start] || ![string is integer $end]} {error "list boundaries must be integers"}
	for {set p $start} {$p <= $end} {incr p} {
	    set d [eval "part $p pr $info"]
	    if {$d != "na"} {puts $channel "\t{$d}"}
	}
    }
    puts $channel "\}"
}

proc blockfile_read_auto_particles {channel read auto} {
    set info "[blockfile $channel read start] [blockfile $channel read toend]"
    set idx 1
    foreach i $info {
	case $i {
	    "id"    { set id $idx; incr idx }
	    "pos"   { set pos $idx; incr idx 3 }
	    "type"  { set type $idx; incr idx }
	    "q"     { set q $idx; incr idx }
	    "f"     { set f $idx; incr idx 3 }
	    "v"     { set v $idx; incr idx 3 }
	    default { error " $i is not a particle property" }
	}
    }
    if {![info exists id] || ![info exists pos]} { error "the fields id and pos are mandatory" }
    set cmd "part \[lindex \$line $id\] \
             pos \[lindex \$line $pos \] \[lindex \$line [expr $pos + 1]\] \[lindex \$line [expr $pos + 2]\]"
    if {[info exists q]} { set cmd "$cmd \
             q \[lindex \$line $q\]" }
    if {[info exists type]} { set cmd "$cmd \
             type \[lindex \$line $type\]" }
    if {[info exists v]} { set cmd "$cmd \
             v  \[lindex \$line $v\] \[lindex \$line [expr $v + 1]\] \[lindex \$line [expr $v + 2]\]"
    }
    if {[info exists f]} { set cmd "$cmd \
             f  \[lindex \$line $f\] \[lindex \$line [expr $f + 1]\] \[lindex \$line [expr $f + 2]\]"
    }
    while {1} {
	set line [blockfile $channel read auto]
	if {[lindex $line 0] != "usertag"} {
	    if {$line != "illstring \}"} { error "particle block ill formed (\"[lindex $line 1]\" unexpected)" }
	    break
	}
	eval $cmd
    }
}

######################################
# bonds support
######################################

proc blockfile_write_bonds {channel write bonds {range "0-end"}} {
    blockfile $channel write start "bonds"

    puts $channel " "
    foreach list $range {
	set list [split $list "-"]
	set len [llength $list]
	if {$len == 1} {
	    if {$list == "all"} {
		set start 0
		set end [setmd max_part]
	    } {
		if {$list == "end"} {set list [setmd max_part]}
		set start $list
		set end   $list
	    }
	} {
	    if {$len != 2} {
		error "list must have form \"<a>-<b> <c>-end...\""
		exit
	    }
	    set start [lindex $list 0]
	    if {$start == "end"} {set start [setmd max_part]}
	    set end [lindex $list 1]
	    if {$end == "end"} {set end [setmd max_part]}
	}
	if {![string is integer $start] || ![string is integer $end]} {error "list boundaries must be integers"}
	for {set p $start} {$p <= $end} {incr p} {
	    set d [eval "part $p pr bonds"]
	    if {$d != "na" && $d != "{}"} {puts $channel "\t{$p $d}"}
	}
    }
    puts $channel "\}"
}

proc blockfile_read_auto_bonds {channel read auto} {
    while {1} {
	set line [blockfile $channel read auto]
	if {[lindex $line 0] != "usertag"} {
	    if {$line != "illstring \}"} { error "particle block ill formed (\"[lindex $line 1]\" unexpected)" }
	    break
	}
	set pid [lindex $line 1] 
	set bl [lindex $line 2]
	foreach inter $bl { part $pid bond delete; eval [concat {part $pid bond} $inter] }
    }
}

######################################
# interactions support
######################################

proc blockfile_write_interactions {channel write interactions} {
    blockfile $channel write start interactions
    puts $channel "\n\t{[join [inter] "\}\n\t\{"]}\n\}"
}

proc blockfile_read_auto_interactions {channel read auto} {
    set data [blockfile $channel read toend]
    foreach d $data { eval "inter $d" }
}

######################################
# random/seed support
######################################

proc blockfile_write_random {channel write random} {
    blockfile $channel write start random
    puts $channel "\n\t{[join [t_random stat] "\} \{"]}\n\}"
}

proc blockfile_read_auto_random {channel read auto} {
    set data [blockfile $channel read toend]
    eval t_random stat [eval concat $data]
}

proc blockfile_write_seed {channel write seed} {
    blockfile $channel write start seed
    puts $channel "\n\t{[join [t_random seed] "\} \{"]}\n\}"
}

proc blockfile_read_auto_seed {channel read auto} {
    set data [blockfile $channel read toend]
    eval t_random seed $data
}

######################################
# setmd variables support
######################################

proc blockfile_write_variable {channel write variable {which "all"}} {
    blockfile $channel write start variable
    if { $which == "all" } {
	puts $channel "\n\t{[join [setmd] "\}\n\t\{"]}"
    } {
	if {[llength $which] == 1} {
	    puts -nonewline $channel " {$which [setmd $which]} "
	} {
	    puts $channel ""
	    foreach wh $which { puts $channel "\t{$wh [setmd $wh]}" }
	}
    }
    puts $channel "\}"
}

# this is more tricky since we want to read old files
proc blockfile_read_auto_variable {channel read auto} {
    set vars [blockfile $channel read toend]
    set vname [lindex $vars 0]
    if {[llength $vname] == 1} {
	# old format
	set data [lindex $vars 1]
	set type [lindex $data 0]
	if { $type != "_ival_" && $type != "_dval_"} { error "old style variable block corrupt"}
	eval "setmd $vname [lrange $data 1 end]"
    } {
	# new format
	foreach vblock $vars {
	    set vname [lindex $vblock 0]
	    set data [lrange $vblock 1 end]
	    if {[catch {eval "setmd $vname $data"} error]} {
		if { $error != "variable is readonly" } {
		    error "blockfile_read_auto_variable: setmd $vname $data reported: $error"
		}
	    }   
	}
    }
}

######################################
# set variables support
######################################

proc blockfile_write_tclvariable {channel write tclvariable {which "all"}} {
    blockfile $channel write start tclvariable
    if { $which == "reallyall" } { set which [info globals] }
    if { $which == "all" } {
	# filter tcl system variables
	set which [info globals]
	foreach i {tcl_version argv argv0 tcl_interactive auto_oldpath errorCode auto_path errorInfo tcl_pkgPath
	    tcl_patchLevel argc tcl_libPath tcl_library} {
	    set ind [lsearch -exact $which $i]
	    set which [lreplace $which $ind $ind]
	}
    }
    if {[llength $which] == 1} {
	global $which
	puts -nonewline $channel " {$which [set $which]} "
    } {
	puts $channel ""
	foreach wh $which { global $wh; if { ! [array exists $wh] } { puts $channel "\t{$wh {[set $wh]}}" } }
    }
    puts $channel "\}"
}

# this is more tricky since we want to read old files
proc blockfile_read_auto_tclvariable {channel read auto} {
    set vars [blockfile $channel read toend]
    puts "--$vars"
    foreach vblock $vars {
	set vname [lindex $vblock 0]
	set data [lrange $vblock 1 end]
	puts "----$vname-$data-"
	global $vname
	if {[catch {eval "set $vname $data"} error]} {
	    if { $error != "" } {
		error "blockfile_read_auto_tclvariable: set $vname $data reported: $error"
	    }
	}   
    }
}

