#
# This set of routines are used to set constraints using the expresso
# constraint feature
#

namespace eval ::mbtools::system_generation {}

proc ::mbtools::system_generation::create_sphericalconstraint { topo box_l args } { 
    mmsg::send [namespace current] "making a spherical constraint "
    # ---- Process Command Line Args --------------------------- #
    
    set options {
	{center.arg  { 0.0 0.0 0.0 } "location of the sphere relative to the box center" }
	{radius.arg  1.0      "radius of the constraint" }
	{dir.arg    1  "constrain out or in -1 is inward" }
	{type.arg  0  "the type of the constraint for the purposes of setting interactions"}

    }
    set usage "Usage: create_uniform_sphere alipid: "
    array set params [::cmdline::getoptions args $options $usage]

    #Setup the center so it is offset from the box center
    set center [list 0.0 0.0 0.0 ]

    #Construct the real center
    lset center 0 [expr [lindex $params(center) 0] + [lindex $box_l 0]/(2.0)]
    lset center 1 [expr [lindex $params(center) 1] + [lindex $box_l 1]/(2.0)]
    lset center 2 [expr [lindex $params(center) 2] + [lindex $box_l 2]/(2.0)]


    constraint sphere center [lindex $center 0] [lindex $center 1] [lindex $center 2] radius $params(radius) direction $params(dir) type $params(type)


    return
    
    


}


