# mathutils
#
#  basic mathematics proceedures that are not included in the rather
#  limited tcllib math package
#

package provide mathutils 0.1.0

namespace eval mathutils {
    namespace export calc_com
    namespace export average
    namespace export stdev
    namespace export acorr
    namespace export normalize
    namespace export scalevec
    namespace export distance
}

source [file join [file dirname [info script]] statistics.tcl]

# ::mathutils::calc_com --
#
# calculate center of mass of a molecule
#
# 
# Arguments:
# m : The molecule id
# topo : the topology
#
proc ::mathutils::calc_com { m  topo } {
    set com " 0.0 0.0 0.0"
    for { set p 1 } { $p < [llength [lindex $topo $m] ] } { incr p } { 
	set pid [lindex [lindex $topo $m] $p]
	set pos [part $pid p p]
	for { set i 0 } { $i < 3 } { incr i } {
	    lset com $i [expr ([lindex $pos $i])/3.0 + [lindex $com $i]]
	}
    }
    return $com
}

# ::mathutils::average --
#
# calculate the average of a list of numerical data
#
proc ::mathutils::average { data { from 0 } { to 0 } } {
	set sum 0
	set cnt 0
	if { $to == 0 || $to >= [llength $data] } {
	    set to [expr [llength $data] - 1 ]
	}
	if { $from < 0 || $from >= [llength $data] } {
	    return 0
	}

	for { set i $from } { $i < [expr $to + 1] } { incr i } {
	    set sum [expr $sum + [lindex $data $i] ]
	    incr cnt
	}


	return [expr $sum/(1.0*$cnt)]
}

# ::mathutils::stdev --
#
# calculate the stdev of a list of numerical data
#
proc ::mathutils::stdev { data { from 0 } { to 0 } } {
	set sum 0
	set cnt 0
	if { $to == 0 || $to >= [llength $data] } {
	    set to [expr [llength $data] - 1 ]
	}
	if { $from < 0 || $from >= [llength $data] } {
	    return 0
	}

	set av [average $data $from $to]

	for { set i $from } { $i < [expr $to + 1] } { incr i } {
	    set sum [expr $sum + ([lindex $data $i] - $av)* ([lindex $data $i] - $av)]
	    incr cnt
	}


	set var [expr $sum/(1.0*$cnt)]
	return [expr sqrt($var)]
}


# ::mathutils::acorr --
#
# calculate an autocorrelation function from a list of numbers
#
proc ::mathutils::acorr { data } {
	set tot [llength $data]
	set to [expr $tot - 1 ]


	set av [average $data ]
	puts $av
	set sum 0
	for { set x 0 } { $x < $tot } { incr x } {
	    set sum [ expr [lindex $data $x]*[lindex $data $x] + $sum]
	}
	set variance [expr $sum/(1.0*$tot) - $av*$av]
	puts $variance

	for { set y 0 } { $y < $tot } { incr y } {
	    set sum1 0
	    set sum2 0
	    set sum3 0
	    for { set x 0 } { $x < [expr $tot - $y] } { incr x } {
		set sum1 [expr [lindex $data $x]*[lindex $data [expr $x + $y]] + $sum1]
		set sum2 [expr $sum2 + [lindex $data $x]]
		set sum3 [expr $sum3 + [lindex $data [expr $x + $y]]]
	    }
#	    puts "$y $sum1 $sum2 $sum3 [expr ($sum1 - $sum2*$sum3)]"
	    set ct [expr ($sum1 - $sum2*$sum3/(1.0*($tot-$y)))/(1.0*($tot - $y)*$variance)]
	    lappend acorr $ct
	}
	return $acorr
}


# ::mathutils::distance --
#
# Distance between two vectors
#
proc ::mathutils::distance { vec1 vec2 } {


    if { [ llength $vec1 ] != [llength $vec2] } {
	mmsg::err [namespace current] "cannot find distance between vectors $vec1 and $vec2 because they are different lengths"
    }

    set sum 0
    for {set i 0 } { $i < [llength $vec1] } { incr i } {
	set diff [expr [lindex $vec1 $i] - [lindex $vec2 $i]]
	set sum [expr $sum + $diff*$diff ]
    }

    return [expr sqrt($sum)]
}


# ::mathutils::normalize --
#
# Normalize a vector
#
proc ::mathutils::normalize { vec } {
    # Find the magnitude
    set mag2 0.0
    foreach val $vec {
	set mag2 [expr $mag2 + $val*$val]
    }
    set mag [expr sqrt($mag2)]
    
    # Divide by the magnitude to produce normalized vector
    foreach val $vec {
	lappend nvec [expr $val/$mag]
    }
    return $nvec
}

# ::mathutils::scalevec --
#
# Scale a vector
#
proc ::mathutils::scalevec { vec scale } {
    foreach val $vec {
	lappend svec [expr $val*$scale]
    }
    return $svec
}


