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
		set end [setmd maxpart]
	    } {
		if {$list == "end"} {set list [setmd maxpart]}
		set start $list
		set end   $list
	    }
	} {
	    if {$len != 2} {
		error "list must have form \"<a>-<b> <c>-end...\""
		exit
	    }
	    set start [lindex $list 0]
	    if {$start == "end"} {set start [setmd maxpart]}
	    set end [lindex $list 1]
	    if {$end == "end"} {set end [setmd maxpart]}
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

proc blockfile_read_particles {channel read particles} {
    set tag [blockfile $channel read start]
    if {$tag != "particles"} { error "blockfile read particles did not find particle block in file $channel" }
    blockfile_read_auto_particles $channel read particles
}
 
proc blockfile_write_bonds {channel write bonds {range "0-end"}} {
    blockfile $channel write start "bonds"

    puts $channel " "
    foreach list $range {
	set list [split $list "-"]
	set len [llength $list]
	if {$len == 1} {
	    if {$list == "all"} {
		set start 0
		set end [setmd maxpart]
	    } {
		if {$list == "end"} {set list [setmd maxpart]}
		set start $list
		set end   $list
	    }
	} {
	    if {$len != 2} {
		error "list must have form \"<a>-<b> <c>-end...\""
		exit
	    }
	    set start [lindex $list 0]
	    if {$start == "end"} {set start [setmd maxpart]}
	    set end [lindex $list 1]
	    if {$end == "end"} {set end [setmd maxpart]}
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

proc blockfile_read_bonds {channel read bonds} {
    set tag [blockfile $channel read start]
    if {$tag != "particles"} { error "blockfile read bonds did not find particle block in file $channel" }
    blockfile_read_auto_bonds $channel read bonds
}
 
proc blockfile_write_interactions {channel write interactions} {
    blockfile $channel write start interactions
    puts $channel "\n\t{[join [inter] "\}\n\t\{"]}\n\}"
}

proc blockfile_read_auto_interactions {channel read auto} {
    set data [blockfile $channel read toend]
    foreach d $data { eval "inter $d" }
}

proc blockfile_read_interactions {channel read interactions} {
    set tag [blockfile $channel read start]
    if {$tag != "particles"} { error "blockfile read interactions did not find interactions block in file $channel" }
    blockfile_read_auto_particles $channel read interactions
}
