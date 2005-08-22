# ::std_analysis::analyze_stray --
#
# Calculate the number of stray lipids
#
namespace eval ::std_analysis {}

proc ::std_analysis::analyze_stray { printflag } {
    variable this
    mmsg::send $this "analyzing number of stray lipids "
    variable av_stray
    variable av_stray_i
    set l_orients [ analyze get_lipid_orients]
#    puts $l_orients
    for { set i 0 } { $i < [llength [lindex $l_orients 1] ] } {incr i } { 
	if { [lindex $l_orients 1 $i ] == 3 || [lindex $l_orients 1 $i ] == 2  } {
	    set av_stray [expr $av_stray + 1]
	}
    }
#    puts "$av_stray $av_stray_i"
    incr av_stray_i

    if { $printflag } {
	mmsg::send $this "strays: [expr $av_stray/(1.0*$av_stray_i)] : $av_stray_i"
	flush stdout
    }
    
}

