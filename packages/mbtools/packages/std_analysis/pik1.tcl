# ::std_analysis::analyze_pik1 --
#
#  Calculate the pressure tensor of the system.  
#

namespace eval ::std_analysis {}

proc ::std_analysis::analyze_pik1 { printflag } {
    variable this
    variable all_particles
    variable av_pik1
    variable av_pik1_i
    mmsg::send $this "analyzing pik1"
    
    set blen [setmd box_l]
    set tot_volume [expr [lindex $blen 0]*[lindex $blen 1]*[lindex $blen 2]]
    
    set ik1_tmp [analyze p_IK1 $tot_volume $all_particles 0]
    
    
    
    set tp 0
    for { set i 0 } { $i < 9 } { incr i } {
	set component  "[lindex [lindex $ik1_tmp $tp] [expr $i+1] ]"
	lset av_pik1  $i [expr [lindex $av_pik1 $i] + $component]
    }
    
    incr av_pik1_i
    
    if { $printflag } {
	set avpx [expr [lindex $av_pik1 0]/($av_pik1_i*1.0)]
	set avpy [expr [lindex $av_pik1 4]/($av_pik1_i*1.0)]
	set avpz [expr [lindex $av_pik1 8]/($av_pik1_i*1.0)]
	mmsg::send $this "<p> : xx:$avpx yy:$avpy yy:$avpz"
	
	set matrixout [format "\n%.5f %.5f %.5f \n%.5f %.5f %.5f \n%.5f %.5f %.5f \n" [lindex $av_pik1 0] [lindex $av_pik1 1] [lindex $av_pik1 2] [lindex $av_pik1 3] [lindex $av_pik1 4] [lindex $av_pik1 5]  [lindex $av_pik1 6] [lindex $av_pik1 7] [lindex $av_pik1 8]]
	mmsg::send $this $matrixout
    }
    mmsg::debug $this "done"
    
}
