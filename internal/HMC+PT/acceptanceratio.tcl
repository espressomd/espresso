#!/bin/sh
source ./hmc-initialize.tcl
source ./hmc-setup.tcl

set start [lindex $argv 1]
set end   [lindex $argv 2]
set step  [lindex $argv 3]

for {set index 39} { $index < $ensemble_num} {incr index} {
    set bjerrum [lindex $bjerrum_list $index]
    puts "BJERRUM $bjerrum"
    exit 0
    if { [file exists "lb_[format %06.3f $bjerrum]/acceptance.txt"] } {
	set lastinput [exec tail -1  "lb_[format %06.3f $bjerrum]/acceptance.txt"]

	set f [open "lb_[format %06.3f $bjerrum]/acceptance.txt" a]
		#puts "CONFIG [format %5d $config] SUCCESS_SWP [format %5d $cum_s_swp] SUCCESS_HMC [format %5d $cum_s_hmc]  SUCC [format %5d $succ]"		
	set cum_s_hmc [lindex $lastinput 5]
	set cum_s_swp [lindex $lastinput 3]
	set cum_f_swp [expr [lindex $lastinput 1] / 10 - $cum_s_swp]
	set cum_f_hmc [expr [lindex $lastinput 1] - $cum_s_hmc]
	puts "[lindex $lastinput 1] $cum_s_hmc $cum_s_swp $cum_f_hmc $cum_f_swp"
    } else {
	set f [open "lb_[format %06.3f $bjerrum]/acceptance.txt" w]
	set cum_f_hmc 0
	set cum_s_hmc 0
	set cum_s_swp 0
	set cum_f_swp 0    
    }

    set fail 0
    set succ 0

    #########################################################################
    set config $start
    set file [open "|gzip -cd lb_[format %06.3f $bjerrum]/$name.t_[format %05d $config].gz" "r"]
    while { [blockfile $file read auto] != "eof" } {}
    if {[catch {close $file} err]} {
	if {![regexp "decompression OK, trailing garbage ignored" $err]} {
	    error $err
	}
    }
    set partold [part 0 print pos]

    for { set config [expr $start+1 ] } { $config <= $end } { incr config $step} {
	set file [open "|gzip -cd lb_[format %06.3f $bjerrum]/$name.t_[format %05d $config].gz" "r"]
	while { [blockfile $file read auto] != "eof" } {}
	if {[catch {close $file} err]} {
	    if {![regexp "decompression OK, trailing garbage ignored" $err]} {
		error $err
	    }
	}
	set partnew [part 0 print pos]
	if { ![expr $config % 10] } {
	    set cum_f_hmc [expr $cum_f_hmc + $fail]
	    set cum_s_hmc [expr $cum_s_hmc + $succ]
	    #puts "HE $fail $succ "
	    if { ![string compare $partold $partnew] } {
		incr cum_f_swp
		#puts "fail $config $partnew $partold "
	    } else {
		incr cum_s_swp		
		#puts "success $config  $partnew $partold" 
	    }
	    if {[lindex $argv 4] == "file" } {
		puts $f "CONFIG [format %5d $config] SUCCESS_SWP [format %5d $cum_s_swp] SUCCESS_HMC [format %5d $cum_s_hmc]  SUCC [format %5d $succ]"		
	    } else {
		puts "CONFIG [format %5d $config] SUCCESS_SWP [format %5d $cum_s_swp] SUCCESS_HMC [format %5d $cum_s_hmc]  SUCC [format %5d $succ]"		
	    }
	    
	    set fail 0
	    set succ 0

	} else {
	    if { ![string compare $partold $partnew] } {
		incr fail
		#puts "fail $config $partnew $partold "
	    } else {
		incr succ		
		#puts "success $config  $partnew $partold" 
	    }
	    #puts "$fail $succ "
	    
	}
	set partold $partnew
    }
    close $f
}
#exit 0
