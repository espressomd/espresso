#  This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
#  It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
#  and by which you are legally bound while utilizing this file in any form or way.
#  There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  You should have received a copy of that license along with this program;
#  if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
#  write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
#  Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#  
######################################
# particle support
######################################

proc blockfile_write_particles {channel write particles {info "id pos v q"} {range "0-end"}} {
    blockfile $channel write start "particles"
    if {![regexp "pos" $info]} {set info "pos $info"}
    if {![regexp "id" $info]} {set info "id $info"}
    if {[regexp "bond" $info]} { error "bonding information cannot be written" }

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
	# the switch checks which property it MIGHT be at all, the second checks wether it is really a prefix
	switch -regexp $i {
	    "^i"      { if {![regexp "^$i" "identity"]} { error " $i is not a particle property" }
		set id $idx; incr idx }
	    "^p"      { if {![regexp "^$i" "position"]} { error " $i is not a particle property" }
		 set pos $idx; incr idx 3 }
	    "^ty"     { if {![regexp "^$i" "type"]} { error " $i is not a particle property" }
		 set type $idx; incr idx }
	    "^m"      { if {![regexp "^$i" "molecule_id"]} { error " $i is not a particle property" }
		 set mol $idx; incr idx }
	    "^q$"     { set q $idx; incr idx }
	    "^fo"     { if {![regexp "^$i" "force"]} { error " $i is not a particle property" }
		 set f $idx; incr idx 3 }
	    "^v"      { if {![regexp "^$i" "v"]} { error " $i is not a particle property" }
		 set v $idx; incr idx 3 }
	    "^qu"     { if {![regexp "^$i" "quat"]} { error " $i is not a particle property" }
		 set quat $idx; incr idx 4 }
	    "^o"      { if {![regexp "^$i" "omega"]} { error " $i is not a particle property" }
		 set omega $idx; incr idx 3 }
	    "^to"     { if {![regexp "^$i" "torque"]} { error " $i is not a particle property" }
		 set torque $idx; incr idx 3 }
	    "^fix"    { if {![regexp "^$i" "fix"]} { error " $i is not a particle property" }
		 set fix $idx; incr idx 3 }
	    "^e"      { if {![regexp "^$i" "ext_force"]} { error " $i is not a particle property" }
		 set ext_force $idx; incr idx 3 }
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
    if {[info exists mol]} { set cmd "$cmd \
             molecule \[lindex \$line $mol\]" }
    if {[info exists v]} { set cmd "$cmd \
             v  \[lindex \$line $v\] \[lindex \$line [expr $v + 1]\] \[lindex \$line [expr $v + 2]\]"
    }
    if {[info exists f]} { set cmd "$cmd \
             f  \[lindex \$line $f\] \[lindex \$line [expr $f + 1]\] \[lindex \$line [expr $f + 2]\]"
    }
    if {[info exists quat]} { set cmd "$cmd \
             quat \[lindex \$line $quat\] \[lindex \$line [expr $quat + 1]\] \[lindex \$line [expr $quat + 2]\] \[lindex \$line [expr $quat + 3]\]"
    }
    if {[info exists omega]} { set cmd "$cmd \
             omega \[lindex \$line $omega\] \[lindex \$line [expr $omega + 1]\] \[lindex \$line [expr $omega + 2]\]"
    }
    if {[info exists torque]} { set cmd "$cmd \
             torque \[lindex \$line $torque\] \[lindex \$line [expr $torque + 1]\] \[lindex \$line [expr $torque + 2]\]"
    }    
    if {[info exists fix]} { set cmd "$cmd \
             fix \[lindex \$line $fix\] \[lindex \$line [expr $fix + 1]\] \[lindex \$line [expr $fix + 2]\]"
    }
    if {[info exists ext_force]} { set cmd "$cmd \
             ext_force \[lindex \$line $ext_force\] \[lindex \$line [expr $ext_force + 1]\] \[lindex \$line [expr $ext_force + 2]\]"
    }    
    while {1} {
	set line [blockfile $channel read auto]
	if {[lindex $line 0] != "usertag"} {
	    if {$line != "illstring \}"} { error "particle block ill formed (\"[lindex $line 1]\" unexpected)" }
	    break
	}
	eval $cmd
    }

    return "particles"
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
	part $pid bond delete
	foreach inter $bl { eval [concat {part $pid bond} $inter] }
    }

    return "bonds"
}

######################################
# interactions support
######################################

proc blockfile_write_interactions {channel write interactions} {
    blockfile $channel write start interactions
    set data [join [inter] "\}\n\t\{"]
    puts $channel "\n\t{$data}\n\}"
}

proc blockfile_read_auto_interactions {channel read auto} {
    set data [blockfile $channel read toend]
    foreach d $data { eval "inter $d" }
}

######################################
# topology support
######################################

proc blockfile_write_topology {channel write topology} {
    blockfile $channel write start topology
    set data [join [analyze set] "\}\n\t\{"]
    puts $channel "\n\t{$data}\n\}"
}

proc blockfile_read_auto_topology {channel read auto} {
    set data [string map {"\n" " "} [blockfile $channel read toend]]
    eval "analyze set $data"
}

######################################
# random/seed support
######################################

# t_random
proc blockfile_write_random {channel write random} {
    blockfile $channel write start random
    puts $channel "\n\t{[join [t_random stat] "\} \{"]}\n\}"
}

proc blockfile_read_auto_random {channel read auto} {
    set data [blockfile $channel read toend]
    eval t_random stat [eval concat $data]

    return "random"
}

proc blockfile_write_seed {channel write seed} {
    blockfile $channel write start seed
    puts $channel "\n\t{[join [t_random seed] "\} \{"]}\n\}"
}

proc blockfile_read_auto_seed {channel read auto} {
    set data [blockfile $channel read toend]
    eval t_random seed $data

    return "seed"
}

# bit_random
proc blockfile_write_bitrandom {channel write bitrandom} {
    blockfile $channel write start bitrandom
    puts $channel "\n\t{[join [bit_random stat] "\} \{"]}\n\}"
}

proc blockfile_read_auto_bitrandom {channel read auto} {
    set data [blockfile $channel read toend]
    eval bit_random stat [eval concat $data]

    return "bitrandom"
}

proc blockfile_write_bitseed {channel write bitseed} {
    blockfile $channel write start bitseed
    puts $channel "\n\t{[join [bit_random seed] "\} \{"]}\n\}"
}

proc blockfile_read_auto_bitseed {channel read auto} {
    set data [blockfile $channel read toend]
    eval bit_random seed $data

    return "bitseed"
}

######################################
# configs support
######################################

proc blockfile_write_configs {channel write configs {range "all"} } {
    blockfile $channel write start configs
    if { "$range" == "all" } { set range 0 } else { set range [expr [analyze stored]-$range] }
    puts -nonewline $channel "\n\t\{ "
    if { ($range >= 0) && ([analyze stored] > 0) } {
	for { set i $range } { $i < [expr [analyze stored]-1] } { incr i } { 
	    puts -nonewline $channel "[analyze configs $i] \}\n\t\{ "
	}
	puts $channel "[analyze configs [expr [analyze stored]-1]] \}\n\}"
    } else { puts $channel "\}\n\}" }
    # puts $channel "\n\t\{[join [analyze configs $i] "\}\n\t\{"]\}\n\}"
}

proc blockfile_read_auto_configs {channel read auto} {
    set data [blockfile $channel read toend]
    foreach d $data { eval "analyze configs $d" }

    return "configs"
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

    return "variable"
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
	foreach ixx {tcl_version argv argv0 tcl_interactive auto_oldpath errorCode 
	    auto_path errorInfo tcl_pkgPath tcl_patchLevel argc tcl_libPath 
	    tcl_library tk_strictMotif tk_library tk_version tk_patchLevel} {
	    set ind [lsearch -exact $which $ixx]
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
#    puts "--$vars"
    foreach vblock $vars {
	set vname [lindex $vblock 0]
	set data [lrange $vblock 1 end]
#	puts "----$vname-$data-"
	global $vname
	if {[catch {eval "set $vname \"$data\""} error]} {
	    if { $error != "" } {
		error "blockfile_read_auto_tclvariable: set $vname $data reported: $error"
	    }
	}   
    }

    return "tclvariable"
}

