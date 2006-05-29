# ::mbtools::analysis::analyze_box_len --
#
# Extract the box dimensions from espresso
# 

namespace eval ::mbtools::analysis {}

namespace eval ::mbtools::analysis::boxl {
    variable av_boxl { 0 0 0 }
    variable av_boxl_i 0
    variable f_tvsbl
    variable verbose

    namespace export setup_boxl
    namespace export analyze_boxl
    namespace export printav_boxl
    namespace export resetav_boxl
}

proc ::mbtools::analysis::boxl::resetav_boxl { } {
    variable av_boxl 
    variable av_boxl_i 
    set av_boxl { 0.0 0.0 0.0 }
    set av_boxl_i 0
}

proc ::mbtools::analysis::boxl::printav_boxl { } {
    variable av_boxl_i
    variable av_boxl
    variable f_tvsbl
    global ::mbtools::analysis::time
    
    if { $av_boxl_i > 0 } {
	set avblx [expr [lindex $av_boxl 0]/($av_boxl_i*1.0)]
	set avbly [expr [lindex $av_boxl 1]/($av_boxl_i*1.0)]
	set avblz [expr [lindex $av_boxl 2]/($av_boxl_i*1.0)]
	puts $f_tvsbl "$time $avblx $avbly $avblz"
    } else {
	::mmsg::warn [namespace current] "can't print average box length"
	flush stdout
    }
    flush $f_tvsbl
}

proc ::mbtools::analysis::boxl::setup_boxl { args } {
    global  ::mbtools::analysis::outputdir
    variable f_tvsbl
    global ::mbtools::analysis::iotype
    global ::mbtools::analysis::suffix
    variable verbose
    ::mmsg::debug [namespace current] "opening $outputdir/time_vs_boxl$suffix "
    
     set options {
	{verbose "print out lots of stuff" }
    }
    set usage "Usage: setup_boxl verbose "
    array set params [::cmdline::getoptions args $options $usage]

    set verbose $params(verbose)

    if { [file exists "$outputdir/time_vs_boxl$suffix"] } {
	set newfile 0
    } else { 
	set newfile 1
    }
    set f_tvsbl [open "$outputdir/time_vs_boxl$suffix" $iotype]
    if { $newfile || $iotype == "w"} {
	puts $f_tvsbl "\# Time boxx boxy boxz"
    }
    

}

proc ::mbtools::analysis::boxl::analyze_boxl {  } {
    ::mmsg::send [namespace current] "analyzing box l"
    variable av_boxl_i
    variable av_boxl
    variable verbose

    set inst_boxl [setmd box_l]
    lset av_boxl 0 [expr [lindex $av_boxl 0] + [lindex $inst_boxl 0] ]
    lset av_boxl 1 [expr [lindex $av_boxl 1] + [lindex $inst_boxl 1] ]
    lset av_boxl 2 [expr [lindex $av_boxl 2] + [lindex $inst_boxl 2] ]

    incr av_boxl_i

    if { $verbose } {
	set avblx [expr [lindex $av_boxl 0]/($av_boxl_i*1.0)]
	set avbly [expr [lindex $av_boxl 1]/($av_boxl_i*1.0)]
	set avblz [expr [lindex $av_boxl 2]/($av_boxl_i*1.0)]
	::mmsg::send [namespace current]  "L: [lindex $inst_boxl 0] [lindex $inst_boxl 1] [lindex $inst_boxl 2] :: <L> $avblx $avbly $avblz"
	flush stdout
    }
    ::mmsg::debug [namespace current] "done"
}
