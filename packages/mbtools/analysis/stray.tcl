# ::mbtools::analysis::analyze_stray --
#
# Calculate the number of stray lipids
#
namespace eval ::mbtools::analysis {}

namespace eval ::mbtools::analysis::stray {

    variable av_stray 0
    variable av_stray_i 0
    variable verbose
    variable f_tvsstray

    namespace export setup_stray
    namespace export analyze_stray
    namespace export printav_stray
    namespace export resetav_stray
}

proc ::mbtools::analysis::stray::resetav_stray { } {
    variable av_stray 
    variable av_stray_i
    set av_stray 0
    set av_stray_i 0 
}

proc ::mbtools::analysis::stray::printav_stray { } {
    variable av_stray
    variable av_stray_i
    global ::mbtools::analysis::time
    variable f_tvsstray

    puts $f_tvsstray "$time [expr $av_stray/(1.0*$av_stray_i)]"
    flush $f_tvsstray
}

proc ::mbtools::analysis::stray::setup_stray { args } {
    global ::mbtools::analysis::outputdir
    global ::mbtools::analysis::suffix
    global ::mbtools::analysis::iotype
    variable f_tvsstray
    variable verbose
    set options {
	{verbose "print out lots of stuff" }
    }
    set usage "Usage: setup_stray verbose "
    array set params [::cmdline::getoptions args $options $usage]

    set verbose $params(verbose)

    if { [file exists "$outputdir/time_vs_stray$suffix"] } {
	set newfile 0
    } else { 
	set newfile 1
		}
    ::mmsg::debug [namespace current] "opening $outputdir/time_vs_stray$suffix "
    set f_tvsstray [open "$outputdir/time_vs_stray$suffix" $iotype ]
    if { $newfile || $iotype == "w"} {
		    puts $f_tvsstray "\# The number of stray lipids vs time"
	puts $f_tvsstray "\# Time num_strays"
    }
    
}


proc ::mbtools::analysis::stray::analyze_stray {  } {
    ::mmsg::send [namespace current] "analyzing number of stray lipids "
    variable av_stray
    variable av_stray_i
    variable verbose
    set l_orients [ analyze get_lipid_orients]
#    puts $l_orients
    for { set i 0 } { $i < [llength [lindex $l_orients 1] ] } {incr i } { 
	if { [lindex $l_orients 1 $i ] == 3 || [lindex $l_orients 1 $i ] == 2  } {
	    set av_stray [expr $av_stray + 1]
	}
    }
#    puts "$av_stray $av_stray_i"
    incr av_stray_i

    if { $verbose } {
	::mmsg::send [namespace current] "strays: [expr $av_stray/(1.0*$av_stray_i)] : $av_stray_i"
	flush stdout
    }
    
}

