# ::std_analysis::analyze_oop --
#
# Calculate the orientational order parameter 
# 

namespace eval ::std_analysis {}

proc ::std_analysis::analyze_oop { printflag } {
    variable this
    mmsg::send $this "analyzing lipid orientation order "
    variable av_oop
    variable av_oop_i
    
    set tmp [analyze lipid_orient_order]

    set av_oop [expr $av_oop + $tmp]
    incr av_oop_i
    
    if { $printflag } {
	mmsg::send $this "s: [expr $av_oop/(1.0*$av_oop_i)]"
	flush stdout
    }
    
    mmsg::debug $this "done"
}

