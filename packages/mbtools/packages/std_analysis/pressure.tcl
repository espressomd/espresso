# ::std_analysis::analyze_pressure --
#
#  Calculate the total pressure of the system and break it into components
#

namespace eval ::std_analysis {}

proc ::std_analysis::analyze_pressure { printflag } {
    variable this
    mmsg::send $this "analyzing pressure"
    variable av_pressure
    variable av_pressure_i
    set p_inst [setmd npt_p_inst_av]
    set p_tmp [analyze pressure]
    lset av_pressure 0 [expr [lindex $av_pressure 0] + $p_inst]
    for { set i 0 } { $i < [llength $p_tmp ] } { incr i } {
	set tmp [lindex $p_tmp $i]
	set ntmp [llength $tmp]	
	if { [ regexp "pressure" $tmp ] } {
	    lset av_pressure 1 [expr [lindex $av_pressure 1] + [lindex $tmp 1]]
	}
	if { [ regexp "ideal" $tmp ] } {
	    lset av_pressure 2 [expr [lindex $av_pressure 2] + [lindex $tmp 1]]
	}
	if { [ regexp "FENE" $tmp ] } {
	    lset av_pressure 3 [expr [lindex $av_pressure 3] + [lindex $tmp 2]]
	}
	if { [ regexp "HARMONIC" $tmp ] } {
	    lset av_pressure 4 [expr [lindex $av_pressure 4] + [lindex $tmp 2]]
	}
	if { [ regexp "nonbonded" $tmp ] } {
	    lset av_pressure 5 [expr [lindex $av_pressure 5] + [lindex $tmp 3]]
	}
	
    }
    incr av_pressure_i
    if { $printflag } {
	mmsg::send $this "p_inst [lindex $av_pressure 0] : pressure [lindex $av_pressure 1] : ideal [lindex $av_pressure 2] : FENE [lindex $av_pressure 3] : HARMONIC [lindex $av_pressure 4] : nonbonded [lindex $av_pressure 5]"	    
    }
    mmsg::debug $this "done"
}
