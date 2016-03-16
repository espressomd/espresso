# Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#   Max-Planck-Institute for Polymer Research, Theory Group
#  
# This file is part of ESPResSo.
#   
# ESPResSo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your
# option) any later version.
#  
# ESPResSo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
# Proceedures for determining the vertical stress profile of a bilayer
#

namespace eval ::mbtools::analysis {}

namespace eval ::mbtools::analysis::stress_profile {
    variable av_stresses
    variable av_stresses_i 0
    variable verbose
    variable zbins
    variable zrange
    namespace export printav_stress_profile
    namespace export setup_stress_profile
    namespace export analyze_stress_profile
    namespace export resetav_stress_profile
}

proc ::mbtools::analysis::stress_profile::resetav_stress_profile { } {
    # Do nothing because we want to average this over the entire simulation
    
}

proc ::mbtools::analysis::stress_profile::printav_stress_profile { } {
    variable f_stressprof
    variable av_stresses_i
    variable av_stresses
    variable zrange
    variable zbins
    global ::mbtools::analysis::outputdir
    global ::mbtools::analysis::suffix

    ::mmsg::send [namespace current] "printing average stress profile"

    set binwidth [expr $zrange/(1.0 * $zbins)]

    if {  $av_stresses_i > 0 } {
	
	set f_stressprof [open "$outputdir/av_zstress$suffix" "w" ]
	# Write a header to the file
	puts -nonewline $f_stressprof "\# zheight "
	for {set i 0} { $i < 3 } {incr i} {
	  for {set j 0} {$j < 3} {incr j} {
	    puts -nonewline $f_stressprof "||$i:$j    "
          }
	}
	puts $f_stressprof ""
	
	# Write out the stress profiles to a file
	for { set bin 0 } { $bin < [llength $av_stresses] } { incr bin } {
	    set currbin [lindex $av_stresses $bin]
	    puts -nonewline $f_stressprof "[expr $bin*$binwidth+($binwidth/2.0)] "
	    for { set i 0 } { $i < [llength $currbin] } { incr i } {
                set numb [format "%1.3e" [expr [lindex $currbin $i]/(1.0*$av_stresses_i)]]
		puts -nonewline $f_stressprof "$numb  "
	    }
	    puts $f_stressprof ""
	}		    
	close $f_stressprof
    }
}

proc ::mbtools::analysis::stress_profile::setup_stress_profile { args } {

    global ::mbtools::analysis::outputdir
    global ::mbtools::analysis::suffix
    global ::mbtools::analysis::iotype
    variable zbins
    variable zrange
    variable verbose

    variable av_stresses
    variable av_stresses_i

    set options {
	{verbose "print out lots of stuff" }
	{zbins.arg "100" "Number of bins used for stress profile analysis" }
	{zrange.arg "10.0" "Range over which to calculate stress profile"}
    }
    set usage "Usage: setup_stress_profile verbose:zbins:zrange "
    array set params [::cmdline::getoptions args $options $usage]

    set verbose $params(verbose)
    set zbins $params(zbins)
    set zrange $params(zrange)
    
    ::mmsg::send [namespace current] "setting up stress profile $zbins bins, and $zrange range"

    #Initialize av_stresses
    set thisbinlist 0.0
    unset thisbinlist
    for { set i 0 } { $i < $zbins + 1 } { incr i } {
	for { set j 0 } { $j < 9 } { incr j } {
	    lappend thisbinlist 0.0
	}
	lappend av_stresses $thisbinlist
	unset thisbinlist
	
    }
    mmsg::send [namespace current] "setup stress profile with $zbins bins and $zrange range"

}

# ::mbtools::analysis::analyze_stress_profile --
# 
# Calculates the stress tensor through the bilayer
#
proc ::mbtools::analysis::stress_profile::analyze_stress_profile { } {
    variable av_stresses
    variable av_stresses_i
    variable zrange
    variable zbins
    variable verbose
    variable av_height
    variable p_tensors

    ::mmsg::send [namespace current] "analyzing stress profile"

    set av_height [average_height]
    ::mmsg::send [namespace current] "average bilayer height is $av_height"
    set stresses [analyze local_stress_tensor 1 1 0 0 0 [expr $av_height - $zrange/2.0] 0 0 $zrange  1 1 $zbins]     
    ::mmsg::send [namespace current] "finished calculating pressure tensors"

    # Bin up the data
    for { set bn 1 } { $bn < $zbins+1 } { incr bn } {
	for { set i 0 } { $i < 9 } { incr i } { 
	    lset av_stresses $bn $i [ expr [lindex $stresses $bn 1 $i] + [lindex $av_stresses $bn $i] ]
	}
    }
    incr av_stresses_i

}

proc ::mbtools::analysis::stress_profile::approximate_centre { args } {
    
    global ::mbtools::analysis::n_particles
    variable slices
    variable part_no
    variable pos
    variable i
    variable no_parts_in_slice_x
    variable no_parts_in_slice_y
    variable no_parts_in_slice_z
    variable box

    set box [setmd box_l]

    set options {
	{slices.arg "8" "number of slices to chop box up into in each direction" }
    }
    set usage "Usage: approximate_centres -slices slices "
    array set params [::cmdline::getoptions args $options $usage]

    set slices $params(slices)
    
    for { set i 0 } { $i < $slices } {incr i} {
	lappend no_parts_in_slice_x 0
	lappend no_parts_in_slice_y 0
	lappend no_parts_in_slice_z 0
    }

    for { set part_no 0 } {$part_no < $n_particles} {incr part_no} {
	set pos [part $part_no print folded_position]
	for {set i 0} {$i <3} {incr i} {
	    if {[lindex $pos 0] < 0 || [lindex $pos 0] > [lindex $box 0]} {
		puts "Positions not folded correctly"
	    }
	}
	variable bin_x
	variable bin_y
	variable bin_z
	set bin_x [expr int(floor([lindex $pos 0]*$slices/[lindex $box 0]))]
	set bin_y [expr int(floor([lindex $pos 1]*$slices/[lindex $box 1]))]
	set bin_z [expr int(floor([lindex $pos 2]*$slices/[lindex $box 2]))]
	lset no_parts_in_slice_x $bin_x [expr [lindex $no_parts_in_slice_x $bin_x] + 1]
	lset no_parts_in_slice_y $bin_y [expr [lindex $no_parts_in_slice_y $bin_y] + 1]
	lset no_parts_in_slice_z $bin_z [expr [lindex $no_parts_in_slice_z $bin_z] + 1]
    }
    
    variable no_parts_in_slice
    lappend no_parts_in_slice $no_parts_in_slice_x
    lappend no_parts_in_slice $no_parts_in_slice_y
    lappend no_parts_in_slice $no_parts_in_slice_z

    variable started
    variable finished
    variable last_hit
    
    variable dir
    variable centre
    set centre ""
    for { set dir 0} {$dir < 3 } {incr dir} {
	set started -1
	set finished -1 
	set last_hit -1
	for { set i 0 } { $i < $slices } {incr i} {
	    if {[lindex $no_parts_in_slice $dir $i] > [expr $n_particles/$slices/2.0]} {
		if  { $last_hit == 0 } {
		    if {$started != -1} {
			puts "Two peaks in density - approximate_centre routine fails"
			return 0
		    } else {
			set started $i
		    }
		}
		set last_hit 1
	    } else {
		if { $last_hit == 1 } {
		    if {$finished != -1} {
			puts "Two peaks in density - approximate_centre routine fails"
			return 0
		    } else {
			set finished $i
		    }
		}
		set last_hit 0
	    }
	
	}
	if { $started == -1 || $finished == -1 } {
	    if { $last_hit == 0 } {
		puts "There's a bug here somewhere"
		return 0
	    } else {
		lappend centre -1
	    }
	} else {
	    if { $started > $finished } {
		lappend centre [expr ($started + ($finished + $slices - $started - 1)/2.0 + 0.5) * [lindex $box 0]/$slices]
	    } else {
		lappend centre [expr ($started + ($finished - $started - 1)/2.0 + 0.5) * [lindex $box 0]/$slices]
	    }
	}
    }
    return $centre
}

proc ::mbtools::analysis::stress_profile::average_height { args } {
    # finds the average height of a membrane - only works if simulations consists of one horizontal membrane
    global ::mbtools::analysis::n_particles
    variable dir

    variable box
    set box [setmd box_l]

    variable centre
    set centre [approximate_centre]
    ::mmsg::send [namespace current] "centre is $centre"
    if { [lindex $centre 0] != -1 || [lindex $centre 1] != -1} {
	puts "Error here"
    }

    variable rough_height
    set rough_height [lindex $centre 2]

    variable part_no
    variable height
    variable height_sum
    set height_sum 0
    for { set part_no 0 } {$part_no < $n_particles} {incr part_no} {
	set height [lindex [part $part_no print folded_position] 2]
	set height [expr $height - $rough_height]
	if { $height < [expr [lindex $box 2]/2.0] } {
	    set $height [expr $height + [lindex $box 2]]
	}
	set height_sum [expr $height_sum + $height]
    }
    return [expr $height_sum/$n_particles + $rough_height]
}
    

