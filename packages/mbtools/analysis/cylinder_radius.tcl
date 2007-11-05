# ::mbtools::analysis::analyze_cylinder_radius --
#
# Extracts the average radius of the cylinder
# plus inner leaflet and outer leaflet
# 
# Structure is av_cylinder_radius { $R_middle $R_inner $R_outer }

namespace eval ::mbtools::analysis {}

namespace eval ::mbtools::analysis::cylinder_radius {
    variable av_cylinder_radius { 0 0 0 }
    variable av_cylinder_radius_i 0
    variable f_tvscr
    variable verbose

    namespace export setup_cylinder_radius
    namespace export analyze_cylinder_radius
    namespace export printav_cylinder_radius
    namespace export resetav_cylinder_radius
    namespace export get_inst_radius
}

proc ::mbtools::analysis::cylinder_radius::resetav_cylinder_radius { } {
    variable av_cylinder_radius 
    variable av_cylinder_radius_i 
    set av_cylinder_radius { 0.0 0.0 0.0 }
    set av_cylinder_radius_i 0
}

proc ::mbtools::analysis::cylinder_radius::printav_cylinder_radius { } {
    variable av_cylinder_radius_i
    variable av_cylinder_radius
    variable f_tvscr
    global ::mbtools::analysis::time
    
    if { $av_cylinder_radius_i > 0 } {
	set avcr  [expr [lindex $av_cylinder_radius 0]/($av_cylinder_radius_i*1.0)]
	set avcri [expr [lindex $av_cylinder_radius 1]/($av_cylinder_radius_i*1.0)]
	set avcro [expr [lindex $av_cylinder_radius 2]/($av_cylinder_radius_i*1.0)]
	puts $f_tvscr "$time $avcr $avcri $avcro"
    } else {
	::mmsg::warn [namespace current] "can't print average radii"
	flush stdout
    }
    flush $f_tvscr
}

proc ::mbtools::analysis::cylinder_radius::setup_cylinder_radius { args } {
    global  ::mbtools::analysis::outputdir
    variable f_tvscr
    global ::mbtools::analysis::iotype
    global ::mbtools::analysis::suffix
    variable verbose
    ::mmsg::debug [namespace current] "opening $outputdir/time_vs_cylinder_radius$suffix "
    
     set options {
	{verbose "print out lots of stuff" }
    }
    set usage "Usage: setup_cylinder_radius verbose "
    array set params [::cmdline::getoptions args $options $usage]

    set verbose $params(verbose)

    if { [file exists "$outputdir/time_vs_cylinder_radius$suffix"] } {
	set newfile 0
    } else { 
	set newfile 1
    }
    set f_tvscr [open "$outputdir/time_vs_cylinder_radius$suffix" $iotype]
    if { $newfile || $iotype == "w"} {
	puts $f_tvscr "\# Time R_middle R_inner R_outer"
    }
    

}

proc ::mbtools::analysis::cylinder_radius::analyze_cylinder_radius {  } {
    ::mmsg::send [namespace current] "analyzing cylinder radius"
    variable av_cylinder_radius_i
    variable av_cylinder_radius
    variable verbose

    set inst_cylinder_radius [::mbtools::analysis::cylinder_radius::get_inst_radius]

    lset av_cylinder_radius 0 [expr [lindex $av_cylinder_radius 0] + [lindex $inst_cylinder_radius 0] ]
    lset av_cylinder_radius 1 [expr [lindex $av_cylinder_radius 1] + [lindex $inst_cylinder_radius 1] ]
    lset av_cylinder_radius 2 [expr [lindex $av_cylinder_radius 2] + [lindex $inst_cylinder_radius 2] ]

    incr av_cylinder_radius_i

    if { $verbose } {
	set avcr  [expr [lindex $av_cylinder_radius 0]/($av_cylinder_radius_i*1.0)]
	set avcri [expr [lindex $av_cylinder_radius 1]/($av_cylinder_radius_i*1.0)]
	set avcro [expr [lindex $av_cylinder_radius 2]/($av_cylinder_radius_i*1.0)]
	::mmsg::send [namespace current]  "R: [lindex $inst_cylinder_radius 0] [lindex $inst_cylinder_radius 1] [lindex $inst_cylinder_radius 2] :: <R> $avcr $avcri $avcro"
	flush stdout
    }
    ::mmsg::debug [namespace current] "done"
}


proc ::mbtools::analysis::cylinder_radius::get_inst_radius {  } {
    # Calculate the instantaneous radii of the cylinder.
    
    # We assume the cylinder is oriented along the z axis.
    # Calculate the position of the axis on the xy plane, 
    # Average over all head beads
    
    set numpart 0
    set cxtmp 0.0
    set cytmp 0.0
    set cztmp 0.0

    set topo [analyze set]
    foreach mol $topo {
	# get the particle number
	set pp [part [lindex $mol 1] print pos]
	set cxtmp [expr $cxtmp + [lindex $pp 0]]
	set cytmp [expr $cytmp + [lindex $pp 1]]
	set cztmp [expr $cztmp + [lindex $pp 2]]
	incr numpart
    }
    # lcenter is the center of mass of all head beads
    set center [list [expr $cxtmp/($numpart*1.0)] [expr $cytmp/($numpart*1.0)] [expr $cztmp/($numpart*1.0)]]
    ::mmsg::debug [namespace current] "Center of cylinder: $center"

    set radinner 0.0
    set radouter 0.0
    
    # Now take the scalar product of the normal vector
    # and the lipid vector.
    # Normal vector: points out of the surface of the cylinder
    # For simplicity and CPU cost, we take the origin starting 
    # from the middle of the z axis.
    # lipid vector: orientation is from tail to head
    set numpart_inner 0
    set numpart_outer 0
    foreach mol $topo {
	# Calculate Normal vector
	set pp1 [part [lindex $mol 1] print pos]
	set cxtmp1 [expr [lindex $pp1 0] - [lindex $center 0]]
	set cytmp1 [expr [lindex $pp1 1] - [lindex $center 1]]
	set cztmp1 [expr [lindex $pp1 2] - [lindex $center 2]]
	set normvector "$cxtmp1 $cytmp1 $cztmp1"
	
	# Calculate lipid vector
	set pp2 [part [lindex $mol 2] print pos]
	set cxtmp2 [expr [lindex $pp1 0] - [lindex $pp2 0]]
	set cytmp2 [expr [lindex $pp1 1] - [lindex $pp2 1]]
	set cztmp2 [expr [lindex $pp1 2] - [lindex $pp2 2]]
	set lipidvector "$cxtmp2 $cytmp2 $cztmp2"
		
	# calculate horizontal distance between origin
	# and lipid.
	# We make sure axis is at the same height as the lipid,
	# in order to get the radius.
	set axis [list [expr [lindex $center 0]] [expr [lindex $center 1]] [expr [lindex $pp1 2]]]

	# Take dot product to determine if it belongs to inner
	# or outer leaflet
	set orientation [::mbtools::utils::dot_product $normvector $lipidvector]
	
	if { $orientation > 0 } {
	    # Outer leaflet
	    set radouter [expr $radouter + [::mbtools::utils::distance $axis $pp1]]
	    incr numpart_outer
	} else {
	    # Inner leaflet
	    set radinner [expr $radinner + [::mbtools::utils::distance $axis $pp1]]
	    incr numpart_inner
	}
    }

    # Average over all radii
    set radouter [expr $radouter/($numpart_outer*1.0)]
    set radinner [expr $radinner/($numpart_inner*1.0)]
    set radmiddle [expr ($radouter + $radinner)*0.5]
    ::mmsg::debug [namespace current] "radii (m/i/o): $radmiddle $radinner $radouter"
    
    return "$radmiddle $radinner $radouter"
}
