# ::std_analysis::analyze_box_len --
#
# Extract the box dimensions from espresso
# 

namespace eval ::std_analysis {}

proc ::std_analysis::analyze_box_len { printflag } {
    variable this 
    mmsg::send $this "analyzing box l"
    variable av_boxl_i
    variable av_boxl
    set inst_boxl [setmd box_l]
    lset av_boxl 0 [expr [lindex $av_boxl 0] + [lindex $inst_boxl 0] ]
    lset av_boxl 1 [expr [lindex $av_boxl 1] + [lindex $inst_boxl 1] ]
    lset av_boxl 2 [expr [lindex $av_boxl 2] + [lindex $inst_boxl 2] ]

    incr av_boxl_i

    if { $printflag } {
	set avblx [expr [lindex $av_boxl 0]/($av_boxl_i*1.0)]
	set avbly [expr [lindex $av_boxl 1]/($av_boxl_i*1.0)]
	set avblz [expr [lindex $av_boxl 2]/($av_boxl_i*1.0)]
	mmsg::send $this  "L: [lindex $inst_boxl 0] [lindex $inst_boxl 1] [lindex $inst_boxl 2] :: <L> $avblx $avbly $avblz"
	flush stdout
    }
    mmsg::debug $this "done"
}
