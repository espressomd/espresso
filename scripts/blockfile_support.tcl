# Copyright (C) 2010,2011,2012,2013 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#  Max-Planck-Institute for Polymer Research, Theory Group
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
#  
######################################
# particle support
######################################

proc blockfile_write_particles {channel write particles {info "id pos v type"} {range "all"}} {
    blockfile $channel write start "particles"
    if {![regexp "pos" $info]} {set info "pos $info"}
    if {![regexp "id" $info]} {set info "id $info"}
    if {[regexp "bond" $info]} { error "bonding information cannot be written" }
    if {[regexp "exclusions" $info]} { error "exclusions cannot be written here, use 'blockfile write exclusions'" }

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
	      "^i"        { if {![regexp "^$i" "identity"]} { error " $i is not a particle property" }
		  set id $idx; incr idx }
	      "^p"        { if {![regexp "^$i" "position"]} { error " $i is not a particle property" }
		  set pos $idx; incr idx 3 }
	      "^ty"       { if {![regexp "^$i" "type"]} { error " $i is not a particle property" }
		  set type $idx; incr idx }
	      "^mo"       { if {![regexp "^$i" "molecule_id"]} { error " $i is not a particle property" }
		  set mol $idx; incr idx }
        "^ma"       { if {![regexp "^$i" "mass"]} { error " $i is not a particle property" }
      set mass $idx; incr idx }
        "^vi"       { if {![regexp "^$i" "virtual"]} { error " $i is not a particle property" }
      set virtual $idx; incr idx }
        "^vs"       { if {![regexp "^$i" "vs_relative"]} { error " $i is not a particle property" }
      set vs_relative $idx; incr idx 2 }
	      "^q$"       { set q $idx; incr idx }
	      "^v"        { if {![regexp "^$i" "v"]} { error " $i is not a particle property" }
		  set v $idx; incr idx 3 }
	      "^di"       { if {![regexp "^$i" "dip"]} { error " $i is not a particle property" }
		  set dip $idx; incr idx 3 }
	      "^qu"       { if {![regexp "^$i" "quat"]} { error " $i is not a particle property" }
		  set quat $idx; incr idx 4 }
	      "^omega_l"  { if {![regexp "^$i" "omega_lab"]} { error " $i is not a particle property" }
		  set omega_lab $idx; incr idx 3 }
	      "^omega_b"  { if {![regexp "^$i" "omega_body"]} { error " $i is not a particle property" }
		  set omega_body $idx; incr idx 3 }
	      "^omega"  { if {![regexp "^$i" "omega"]} { error " $i is not a particle property" }
		  set omega $idx; incr idx 3 }
	      "^torque_l" { if {![regexp "^$i" "torque_lab"]} { error " $i is not a particle property" }
		  set torque_lab $idx; incr idx 3 }
	      "^torque_b" { if {![regexp "^$i" "torque_body"]} { error " $i is not a particle property" }
		  set torque_body $idx; incr idx 3 }
	      "^torque" { if {![regexp "^$i" "torque"]} { error " $i is not a particle property" }
		  set torque $idx; incr idx 3 }
	      "^tbf" { if {![regexp "^$i" "tbf"]} { error " $i is not a particle property" }
		  set tbf $idx; incr idx 3 }
	      "^fix"      { if {![regexp "^$i" "fix"]} { error " $i is not a particle property" }
		  set fix $idx; incr idx 3 }
	      "^f"        { if {![regexp "^$i" "force"]} { error " $i is not a particle property" }
		  set f $idx; incr idx 3 }
	      "^ext_f"    { if {![regexp "^$i" "ext_force"]} { error " $i is not a particle property" }
		  set ext_force $idx; incr idx 3 }
	      "^ext_t"    { if {![regexp "^$i" "ext_torque"]} { error " $i is not a particle property" }
		  set ext_torque $idx; incr idx 3 }
	      default { error "$i is not a particle property or it is not supported by the blockfile mechanism" }
	  }
  }
  if {![info exists id] || ![info exists pos]} { error "the fields id and pos are mandatory" }

  set cmd "part \[lindex \$line $id\] \
           pos \[lindex \$line $pos \] \[lindex \$line [expr $pos + 1]\] \[lindex \$line [expr $pos + 2]\]"
  if {[info exists q]} { set cmd "$cmd \
           q \[lindex \$line $q\]" }
  if {[info exists dip]} { set cmd "$cmd \
           dip  \[lindex \$line $dip\] \[lindex \$line [expr $dip + 1]\] \[lindex \$line [expr $dip + 2]\]"
  }
  if {[info exists type]} { set cmd "$cmd \
           type \[lindex \$line $type\]" }
  if {[info exists mol]} { set cmd "$cmd \
           molecule \[lindex \$line $mol\]" }
  if {[info exists mass]} { set cmd "$cmd \
           mass \[lindex \$line $mass\]" }
  if {[info exists virtual]} { set cmd "$cmd \
           virtual \[lindex \$line $virtual\]" }
  if {[info exists vs_relative]} { set cmd "$cmd \
           vs_relative \[lindex \$line $vs_relative\] \[lindex \$line [expr $vs_relative +1]\]" }
  if {[info exists v]} { set cmd "$cmd \
           v  \[lindex \$line $v\] \[lindex \$line [expr $v + 1]\] \[lindex \$line [expr $v + 2]\]"
  }
  if {[info exists f]} { set cmd "$cmd \
           f  \[lindex \$line $f\] \[lindex \$line [expr $f + 1]\] \[lindex \$line [expr $f + 2]\]"
  }
  if {[info exists quat]} { set cmd "$cmd \
           quat \[lindex \$line $quat\] \[lindex \$line [expr $quat + 1]\] \[lindex \$line [expr $quat + 2]\] \[lindex \$line [expr $quat + 3]\]"
  }
  if {[info exists omega_lab]} { set cmd "$cmd \
           omega_lab \[lindex \$line $omega_lab\] \[lindex \$line [expr $omega_lab + 1]\] \[lindex \$line [expr $omega_lab + 2]\]"
  }
  if {[info exists omega_body]} { set cmd "$cmd \
           omega_body \[lindex \$line $omega_body\] \[lindex \$line [expr $omega_body + 1]\] \[lindex \$line [expr $omega_body + 2]\]"
  } 
  if {[info exists omega]} { set cmd "$cmd \
           omega_body \[lindex \$line $omega\] \[lindex \$line [expr $omega + 1]\] \[lindex \$line [expr $omega + 2]\]"
  }   
  if {[info exists torque_lab]} { set cmd "$cmd \
           torque_lab \[lindex \$line $torque_lab\] \[lindex \$line [expr $torque_lab + 1]\] \[lindex \$line [expr $torque_lab + 2]\]"
  }  
  if {[info exists torque_body]} { set cmd "$cmd \
           torque_body \[lindex \$line $torque_body\] \[lindex \$line [expr $torque_body + 1]\] \[lindex \$line [expr $torque_body + 2]\]"
  }  
  if {[info exists torque]} { set cmd "$cmd \
           torque_body \[lindex \$line $torque\] \[lindex \$line [expr $torque + 1]\] \[lindex \$line [expr $torque + 2]\]"
  }  
  if {[info exists tbf]} { set cmd "$cmd \
           torque_body \[lindex \$line $tbf\] \[lindex \$line [expr $tbf + 1]\] \[lindex \$line [expr $tbf + 2]\]"
  }
  if {[info exists fix]} { set cmd "$cmd \
           fix \[lindex \$line $fix\] \[lindex \$line [expr $fix + 1]\] \[lindex \$line [expr $fix + 2]\]"
  }
  if {[info exists ext_force]} { set cmd "$cmd \
           ext_force \[lindex \$line $ext_force\] \[lindex \$line [expr $ext_force + 1]\] \[lindex \$line [expr $ext_force + 2]\]"
  }
  if {[info exists ext_torque]} { set cmd "$cmd \
           ext_torque \[lindex \$line $ext_torque\] \[lindex \$line [expr $ext_torque + 1]\] \[lindex \$line [expr $ext_torque + 2]\]"
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

proc blockfile_write_bonds {channel write bonds {range "all"}} {
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
    return "interactions"
}

######################################
# constraints support
######################################
proc blockfile_write_constraints {channel write constraints} {
    blockfile $channel write start constraints
    set data [join [constraint] "\}\n\t\{"]
    puts $channel "\n\t{$data}\n\}"
}

proc blockfile_read_auto_constraints {channel read auto} {
    set data [blockfile $channel read toend]
    foreach d $data { eval "constraint [lrange $d 1 end]" }
    return "constraints"
}

######################################
# exclusions support
######################################

proc blockfile_write_exclusions {channel write exclusions {range "0-end"}} {
    blockfile $channel write start "exclusions"

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
	    set d [eval "part $p pr exclusions"]
	    if {$d != "na" && $d != "{}"} {puts $channel "\t{$p $d}"}
	}
    }
    puts $channel "\}"
}

proc blockfile_read_auto_exclusions {channel read auto} {
    while {1} {
	set line [blockfile $channel read auto]
	if {[lindex $line 0] != "usertag"} {
	    if {$line != "illstring \}"} { error "particle block ill formed (\"[lindex $line 1]\" unexpected)" }
	    break
	}
	set pid [lindex $line 1]
	set len [expr [llength $line] -1]
	for {set i 2} {$i <= $len} {incr i} {
	    set excid [lindex $line $i]
	    eval [concat { part $pid exclude $excid} ]
	}
    }

    return "exclusions"
}

######################################
# integrate support
######################################

proc blockfile_write_integrate {channel write integrate} {
    blockfile $channel write start integrate
    set data [join [integrate] "\}\n\t\{"]
    puts $channel "\n\t{$data}\n\}"
}

proc blockfile_read_auto_integrate {channel read auto} {
    set data [blockfile $channel read toend]
    foreach d $data { eval "integrate $d" }
    return "integrate"
}

######################################
# thermostat support
######################################

proc blockfile_write_thermostat {channel write thermostat} {
    blockfile $channel write start thermostat
    set data [join [thermostat] "\}\n\t\{"]
    puts $channel "\n\t\{ off \}\n\t{$data}\n\}"
}

proc blockfile_read_auto_thermostat {channel read auto} {
    set data [blockfile $channel read toend]
    foreach d $data { eval "thermostat $d" }
    return "thermostat"
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
    return "topology"
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
	    set which [lindex $which 0]
	    puts -nonewline $channel " {$which [setmd $which]} "
	} {
	    puts $channel ""
	    foreach wh $which { puts $channel "\t{$wh [setmd $wh]}" }
	}
    }
    puts $channel "\}"
}

proc blockfile_read_auto_variable {channel read auto} {
    global blockfile_variable_blacklist blockfile_variable_whitelist
    set vars [blockfile $channel read toend]
    # new format
    foreach vblock $vars {
	set vname [lindex $vblock 0]
	if {[info exists blockfile_variable_blacklist] && \
		[lsearch -exact $blockfile_variable_blacklist $vname] != -1} { continue }
	if {[info exists blockfile_variable_whitelist] && \
		[lsearch -exact $blockfile_variable_whitelist $vname] == -1} { continue }
	set data [lrange $vblock 1 end]
	if {[catch {eval "setmd $vname $data"} error]} {
	    switch -glob $error {
		"min_num_cells must be at least*" {
		    puts stderr "WARNING: min_num_cells incompatible with current n_nodes, ignoring"
		    continue
		}
		"node grid does not fit n_nodes" {
		    puts stderr "WARNING: node_grid incompatible with current n_nodes, ignoring"
		    continue
		}
		"variable is readonly*" { continue }
	    }
	    error "blockfile_read_auto_variable: setmd $vname $data reported: $error"
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
	set which [lindex $which 0]
	global $which
	if { ![array exists $which] } {
	    puts -nonewline $channel " {$which [set $which]} "
	} { puts -nonewline $channel " {array $which [array get $which]} " }
    } {
	puts $channel ""
	foreach wh $which { global $wh; if { ! [array exists $wh] } { puts $channel "\t{$wh [set $wh]}" } }
    }
    puts $channel "\}"
}

proc blockfile_read_auto_tclvariable {channel read auto} {
    global blockfile_tclvariable_blacklist blockfile_tclvariable_whitelist
    set vars [blockfile $channel read toend]
    foreach vblock $vars {
	set vname [lindex $vblock 0]
	if {[info exists blockfile_tclvariable_blacklist] && \
		[lsearch -exact $blockfile_tclvariable_blacklist $vname] != -1} { continue }
	if {[info exists blockfile_tclvariable_whitelist] && \
		[lsearch -exact $blockfile_tclvariable_whitelist $vname] == -1} { continue }
	set data [lrange $vblock 1 end]
	if { "$vname" != "array" } {
	    global $vname
	    if {[catch {eval "set $vname {$data}"} error]} {
		if { $error != "" } {
		    error "blockfile_read_auto_tclvariable: set $vname $data reported: $error"
		}
	    }   
	} else { 
	    set vname [lindex $vblock 1]
	    set data [lrange $vblock 2 end]
	    global $vname
	    if {[catch {eval "array set $vname {$data}"} error]} {
		if { $error != "" } {
		    error "blockfile_read_auto_tclvariable: set $vname $data reported: $error"
		}
	    }   
	}
    }

    return "tclvariable"
}

