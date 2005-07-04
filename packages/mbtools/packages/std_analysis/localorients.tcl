#
# Procedures for determining the distribution of local lipid orientations
#

namespace eval ::std_analysis {}

# ::std_analysis::analyze_localorients --
#  
# Calculates the projection of the lipid orientation vector onto the
# xy plane for each lipid and then bins the absolute values of these
# vectors.
#
proc ::std_analysis::analyze_localorients { printflag } {
    variable this
    variable av_localorients
    variable av_localorients_i
    variable localorientsnbins
    variable lorange
    mmsg::send $this "analyzing local orients"

    set nbins $localorientsnbins
    set binwidth [expr $lorange/($nbins*1.0)]

    # Internal copy of l_orients ( same molecule order as topology )
    set l_orients [lindex [ analyze get_lipid_orients ] 1]

    set localorients [lindex [analyze lipid_orient_order all] 1]

    foreach or $localorients {
	set x [lindex $or 0]
	set y [lindex $or 1]

	# First figure out the absolute value of x and y components together
	set xyproj [expr sqrt($x*$x + $y*$y)]

	# Work out the bin to which this belongs
	if { $xyproj < $lorange } {
	    set thisbin [expr int(floor(($xyproj)/$binwidth))]
	    # Increment the relevant bin
	    lset av_localorients $thisbin [expr [lindex $av_localorients $thisbin] + 1]
	}



    }

    incr av_localorients_i
    
}


