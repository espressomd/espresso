#
# Proceedures for analyzing the lipid flip-flop rate
#

# ::std_analysis::Compliment --
#
# Function called by "analyze_flipflop".  Compliment takes two lists
# of particle orientations and compares them to find the number of
# particles whose orientation is different between the two.
#

namespace eval ::std_analysis {}

proc ::std_analysis::Compliment { list1 list2 } {
    set flip 0
    set bilcnt 0
    for { set i 0 } { $i < [llength [lindex $list1 1] ] } {incr i } { 
	#	    puts "[lindex $list1 1 $i ] : [lindex $list2 1 $i]"
	if { [lindex $list1 1 $i ] != [lindex $list2 1 $i] && [lindex $list1 1 $i ] != 2 && [lindex $list2 1 $i] != 2 && [lindex $list1 1 $i ] != 3 && [lindex $list2 1 $i] != 3 } {
	    incr flip
	    
	    #	puts "flip $flip : [lindex $list1 1 $i ] : [lindex $list2 1 $i] "
	    #	flush stdout
	}
	if { [lindex $list1 1 $i ] != 2 && [lindex $list2 1 $i] != 2 && [lindex $list1 1 $i ] != 3 && [lindex $list2 1 $i] != 3 } {
	    incr bilcnt
	}
    }
    return [expr ($bilcnt - $flip)/(1.0*$bilcnt)]
    
}

# ::std_analysis::analyze_flipflop --
#
#  Obtains the flip state of the current system and compares it to the
#  starting state to determine the number of flips that have occurred.
#
proc ::std_analysis::analyze_flipflop { printflag } {
    variable this
    variable av_flip
    variable av_flip_i
    variable f_tvsflip
    variable l_orients
    variable l_orients_start

    mmsg::send $this "analyzing flip-flop"

    set l_orients [ analyze get_lipid_orients ]
    
    set flip [ Compliment $l_orients $l_orients_start ]
    set av_flip [expr $av_flip + $flip ]
    incr av_flip_i
    
    if { $printflag } { 
	mmsg::send $this "flipflip: $flip <flipflop>: [expr $av_flip/($av_flip_i*1.0) ]"
	flush stdout
    }
}


