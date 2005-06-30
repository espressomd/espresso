#
# Procedures for determining the distribution of local relative lipid heights
#

namespace eval ::std_analysis {}


proc ::std_analysis::foldpart { pos box } {
    for { set i 0 } { $i < 3 } { incr i } {
	set x  [lindex $pos $i]
	set x [expr $x - floor($x/[lindex $box $i])*[lindex $box $i] ]
	lset pos $i $x
    }
    return $pos
}

proc ::std_analysis::nbhoodlipids { l_orients orient headp rcatch molid } {
    # Search for the neighbours around this head bead
    set nblist [analyze nbhood planar 1 1 0 $headp $rcatch  ]

    # Turn the bead neighbourlist into a molecule neighbour list
    set molnblist -1 
    foreach bead $nblist {
	set thismol [part $bead print mol]
	# If the molecule id is not yet in our mol list then add it
	# We also check if the lipid is in the same monolayer and only include it if it is
	
	#		puts "$orient [lindex $l_orients $thismol]"
	if { [lsearch $molnblist $thismol] == -1 && $thismol != $molid && [lindex $l_orients $thismol] == $orient } {
	    lappend molnblist $thismol
	    if { [lindex $molnblist 0 ] == -1 } {
		set molnblist $thismol
	    }		
	}
    }	    
    return $molnblist
}


proc ::std_analysis::calc_localheights {} {
    variable topology
    variable rcatch
    set cmx 0
    set cmy 0
    set cmz 0

    set box [setmd box_l]

    # Internal copy of l_orients ( same molecule order as topology )
    set l_orients [lindex [ analyze get_lipid_orients ] 1]

    # Check that hdiffs is uninitialized
    set hdiffs 0
    unset hdiffs

    set molid 0
    foreach mol $topology {

	set headp [lindex $mol 1]
	set tb [lindex $mol end ]
	# For some old simulations where 0 was a special hb type we need this
	if { [part $headp print type] != 0 && [part $tb print type] == 0 } {
	    set headp $tb
	}
	set orient [lindex $l_orients $molid]

	# Do not perform any analysis for stray beads
	if { $orient < 2 } {

#	    puts "[part $headp print pos]"
	    set hpos [foldpart [part $headp print pos] $box ]

#	    puts "$molid $headp"


	    set tries 0
	    set maxtries 100
	    set molnblist ""
	    while { [llength $molnblist] != 6 && $tries < $maxtries } {
		if { [llength $molnblist] < 6 } {
		    set rcatch [expr $rcatch + 0.01 ]
		} else {
		    set rcatch [expr $rcatch - 0.01 ]
		}

		if { $rcatch <= 0 } {
		    ::mmsg::err [namespace current] "catch radius went to zero when trying to find lipid neighbours"
		}

		set molnblist [nbhoodlipids $l_orients $orient $headp $rcatch $molid ]
#		puts "[llength $molnblist] $tries $rcatch"
		incr tries

	    }

	    # Use the molnblist to work out height differences
	    foreach nbmol $molnblist {
		set hb [lindex [lindex $topology $nbmol] 1 ]
		set tb [lindex [lindex $topology $nbmol] end ]
		# For some old simulations where 0 was a special hb type we need this
		if { [part $hb print type] != 0 && [part $tb print type] == 0 } {
		    set hb $tb
		}
		set hbpos [foldpart [part $hb print pos] $box]
		lappend hdiffs [expr [lindex $hpos 2] - [lindex $hbpos 2]]
#		puts "[lindex $topology $nbmol] [part $hb print type] $hpos $hbpos $hb"
		
	    }
	    
#	    puts "$molnblist $hdiffs"
	    
	}
	incr molid
    }

    return $hdiffs

}

# ::std_analysis::analyze_localheights --
# 
# Calculates the relative heights of all nearest neighbours of each lipid and bins these into a histogram
#
proc ::std_analysis::analyze_localheights { printflag } {
    variable this
    variable av_localheights
    variable av_localheights_i
    variable localheightsnbins
    variable lhrange

    set nbins $localheightsnbins
    set binwidth [expr $lhrange/($nbins*1.0)]

    mmsg::send $this "analyzing local heights"

    set localheights [calc_localheights]

    # Bin up the localheights
    foreach height $localheights {

	set thisbin [expr int(floor(($height+$lhrange/(2.0))/$binwidth))]
#	puts "$height $thisbin"
	if { $thisbin > 0 && $thisbin < $nbins } {
#	    puts "$thisbin [llength $av_localheights]"
	    lset av_localheights $thisbin [expr [lindex $av_localheights $thisbin] + 1 ]
	}
    }

    incr av_localheights_i
    
}


