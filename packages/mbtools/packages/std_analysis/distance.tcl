# ::std_analysis::analyze_box_len --
#
# Extract the box dimensions from espresso
# 

namespace eval ::std_analysis {}

proc ::std_analysis::analyze_distance { printflag } {
    variable this 
    mmsg::send $this "analyzing the distance between 2 particles"
    variable av_dist_i
    variable av_dist

    variable middlebead
    
    set middlebead [::system_generation::get_middlebead ]

    set inst_dist [bond_length [lindex $middlebead 0] [lindex $middlebead 1]] 

    lset av_dist 0 [expr [lindex $av_dist 0] + [lindex $inst_dist 0] ]

    #lset av_dist 1 [expr [lindex $av_dist 1] + [lindex $inst_dist 1] ]
    #lset av_dist 2 [expr [lindex $av_dist 2] + [lindex $inst_dist 2] ]

    incr av_dist_i

  
    if { $printflag } {
	set avdistx [expr [lindex $av_dist 0]/($av_dist_i*1.0)]

	#set avdisty [expr [lindex $av_dist 1]/($av_dist_i*1.0)]
	#set avdistz [expr [lindex $av_dist 2]/($av_dist_i*1.0)]



	mmsg::send $this  "L: [lindex $inst_dist 0] :: <L> $avdistx"
	flush stdout
    }
    mmsg::debug $this "done"
}
