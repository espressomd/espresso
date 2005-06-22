#
# Proceedures for determining the vertical density profile of a bilayer
#

namespace eval ::std_analysis {}

# ::std_analysis::analyze_density_profile --
# 
# Calculates the density of each lipid bead type as a function of
# vertical distance relative to the bilayer midplane
#
proc ::std_analysis::analyze_density_profile { printflag } {
    variable this
    variable av_densities
    variable av_densities_i
    variable hrange
    variable beadtypes
    variable nbins
    variable profilenogrid

    mmsg::send $this "analyzing density profile"

    if { $profilenogrid } {
	set densities [ analyze bilayer_density_profile $hrange $nbins $beadtypes nogrid ]
    } else {
	set densities [ analyze bilayer_density_profile $hrange $nbins $beadtypes  ]
    }

    set densities [lindex $densities 1]

    for { set bn 0 } { $bn < $nbins } { incr bn } {
	for { set bt 0 } { $bt < [expr 2*[llength $beadtypes]] } { incr bt } {
	    lset av_densities $bn $bt [expr [lindex $densities $bt $bn] + [lindex $av_densities $bn $bt]]
	}
    }
    incr av_densities_i
    
}


