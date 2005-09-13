# ::mbtools::analysis::analyze_oop --
#
# Calculate the orientational order parameter 
# 

namespace eval ::mbtools::analysis {}

namespace eval ::mbtools::analysis::orient_order {
    variable av_oop 0.0
    variable av_oop_i 0
    variable f_tvsoop
    variable verbose
    namespace export setup_orient_order
    namespace export analyze_orient_order
    namespace export printav_orient_order
    namespace export resetav_orient_order
}

proc ::mbtools::analysis::orient_order::resetav_orient_order { } {
    variable av_oop 
    variable av_oop_i 
    set av_oop 0
    set av_oop_i 0
}

proc ::mbtools::analysis::orient_order::printav_orient_order { } {
    variable av_oop
    variable av_oop_i
    variable f_tvsoop
    global ::mbtools::analysis::time

    puts $f_tvsoop "$time [expr $av_oop/(1.0*$av_oop_i)]"
    flush $f_tvsoop
}

proc ::mbtools::analysis::orient_order::setup_orient_order { args } {
    variable f_tvsoop
    global ::mbtools::analysis::outputdir
    global ::mbtools::analysis::iotype
    global ::mbtools::analysis::suffix
    variable verbose
    set options {
	{verbose "print out lots of stuff" }
    }
    set usage "Usage: setup_orient_order verbose "
    array set params [::cmdline::getoptions args $options $usage]

    set verbose $params(verbose)


    if { [file exists "$outputdir/time_vs_oop$suffix"] } {
	set newfile 0
    } else { 
	set newfile 1
    }
    ::mmsg::debug [namespace current] "opening $outputdir/time_vs_oop$suffix "
    set f_tvsoop [open "$outputdir/time_vs_oop$suffix" $iotype ]
    if { $newfile || $iotype == "w"} {
	puts $f_tvsoop "\# S = 0.5*<3*(a_i.n)^2 -1>i"
	puts $f_tvsoop "\# where a_i is the orientation of the ith lipid and n is the average bilayer normal"
	puts $f_tvsoop "\# Time S"
    }
}


proc ::mbtools::analysis::orient_order::analyze_orient_order {  } {
    ::mmsg::send [namespace current] "analyzing lipid orientation order "
    variable av_oop
    variable av_oop_i
    variable verbose
    set tmp [analyze lipid_orient_order]

    set av_oop [expr $av_oop + $tmp]
    incr av_oop_i
    
    if { $verbose } {
	::mmsg::send [namespace current] "s: [expr $av_oop/(1.0*$av_oop_i)]"
	flush stdout
    }
    
    ::mmsg::debug [namespace current] "done"
}

