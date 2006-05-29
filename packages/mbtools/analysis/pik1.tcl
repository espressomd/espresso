# ::mbtools::analysis::analyze_pik1 --
#
#  Calculate the pressure tensor of the system.  
#

namespace eval ::mbtools::analysis {}

namespace eval ::mbtools::analysis::pik1 {
    variable all_particles
    variable f_tvspik1
    variable verbose
    variable av_pik1 { 0 0 0 0 0 0 0 0 0 }
    variable av_pik1_i 0
    namespace export setup_pik1
    namespace export analyze_pik1
    namespace export printav_pik1
    namespace export resetav_pik1
}

proc ::mbtools::analysis::pik1::resetav_pik1 { } {
    variable av_pik1
    variable av_pik1_i
    set av_pik1 { 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 }
    set av_pik1_i 0
}

proc ::mbtools::analysis::pik1::printav_pik1 { } {
    global ::mbtools::analysis::time
    variable av_pik1
    variable f_tvspik1
    variable av_pik1_i
    variable verbose

    if { $av_pik1_i > 0 } {
	puts -nonewline $f_tvspik1 "$time "
	for { set v 0 } { $v < [llength $av_pik1] } {incr v} {
	    puts -nonewline $f_tvspik1 "[expr [lindex $av_pik1 $v]/($av_pik1_i*1.0)] "
	    
	}
	puts $f_tvspik1 " "
    } else {
	mmsg::warn [namespace current] "can't print average pik1"
	flush stdout		       
    }
    flush $f_tvspik1
}

proc ::mbtools::analysis::pik1::setup_pik1 { args } {
    # Use the global command to access variables in the std_analysis namespace
    global ::mbtools::analysis::n_particles
    global ::mbtools::analysis::suffix
    global ::mbtools::analysis::outputdir
    global ::mbtools::analysis::iotype

    variable verbose
    variable all_particles
    variable f_tvspik1

    set options {
	{verbose "print out lots of stuff" }
    }
    set usage "Usage: setup_pik1 verbose "
    array set params [::cmdline::getoptions args $options $usage]

    set verbose $params(verbose)

    for { set j 0 } { $j < $n_particles } { incr j } {
	lappend all_particles $j
    }
    mmsg::debug [namespace current] "opening $outputdir/time_vs_pik1$suffix "
    
    if { [file exists "$outputdir/time_vs_pik1$suffix"] } {
	set newfile 0
    } else { 
	set newfile 1
    }
    set f_tvspik1 [open "$outputdir/time_vs_pik1$suffix" $iotype]
    if { $newfile || $iotype == "w"} {
	puts $f_tvspik1 "\# Components of the total pressure tensor in row major order"
	puts $f_tvspik1 "\# pxx pxy pxz pyx pyy pyz pzx pzy pzz"
    }
}

proc ::mbtools::analysis::pik1::analyze_pik1 {  } {
    variable all_particles
    variable av_pik1
    variable av_pik1_i
    variable verbose
    mmsg::send [namespace current] "analyzing pik1"
    
    set blen [setmd box_l]
    set tot_volume [expr [lindex $blen 0]*[lindex $blen 1]*[lindex $blen 2]]
    
    set ik1_tmp [analyze p_IK1 $tot_volume $all_particles 0]
    
    
    
    set tp 0
    for { set i 0 } { $i < 9 } { incr i } {
	set component  "[lindex [lindex $ik1_tmp $tp] [expr $i+1] ]"
	lset av_pik1  $i [expr [lindex $av_pik1 $i] + $component]
    }
    
    incr av_pik1_i
    
    if { $verbose } {
	set avpx [expr [lindex $av_pik1 0]/($av_pik1_i*1.0)]
	set avpy [expr [lindex $av_pik1 4]/($av_pik1_i*1.0)]
	set avpz [expr [lindex $av_pik1 8]/($av_pik1_i*1.0)]
	mmsg::send [namespace current] "<p> : xx:$avpx yy:$avpy yy:$avpz"
	
	set matrixout [format "\n%.5f %.5f %.5f \n%.5f %.5f %.5f \n%.5f %.5f %.5f \n" [lindex $av_pik1 0] [lindex $av_pik1 1] [lindex $av_pik1 2] [lindex $av_pik1 3] [lindex $av_pik1 4] [lindex $av_pik1 5]  [lindex $av_pik1 6] [lindex $av_pik1 7] [lindex $av_pik1 8]]
	mmsg::send [namespace current] $matrixout
    }
    mmsg::debug [namespace current] "done"
    
}
