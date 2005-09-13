
# mbtools::utils
#
#  basic mathematics proceedures that are not included in the rather
#  limited tcllib math package
#

namespace eval mbtools::utils {
    namespace export uniquelist
    namespace export average
    namespace export stdev
    namespace export acorr
    namespace export normalize
    namespace export scalevec
    namespace export distance
    namespace export find_proportions
    namespace export matrix_multiply
}


# ::mbtools::utils::dot_product
proc ::mbtools::utils::dot_product { A B } {
    set dp 0
    for { set i 0 } { $i < 3 } { incr i } {
	set dp [expr $dp + [lindex $A $i]*[lindex $B $i] ]
    }
    return $dp

}



# ::mbtools::utils::matrix_vec_multiply --
# 
# Perform the multiplication of a matrix A and a vector B where A is the
# first argument and B is the second argument.  The routine will
# return AxB
# 
# A should be formatted as a tcl list where each element of the
# list represents the elements in each row of the matrix
#

proc ::mbtools::utils::matrix_vec_multiply { A B } {
    set AB 0
    unset AB

    set side [llength $A]

    for { set i 0 } { $i < $side } { incr i } {
	set sum 0
	# index the row vector from A
	set RA [lindex $A $i]

	# Now find the dot product of this row with B
	for { set j 0 } { $j < $side } { incr j } {
	    set sum [expr $sum + [lindex $RA $j]*[lindex $B $j] ]
	}
	lappend AB $sum
    }
    return $AB
}



# ::mbtools::utils::calc_proportions --
#
# Calculate the proportion of each integer in a list of integers
#
# 
# Arguments:
# ilist: the list of integers
#
proc ::mbtools::utils::calc_proportions { ilist } {

    # Find all the different integers
    set itypelist [lindex $ilist 0]
    lappend plist [list $itypelist 0]
    foreach i $ilist {
	if { [lsearch $itypelist $i] == -1 } {
	    lappend itypelist $i
	    lappend plist [list $i 0 ]
	}
    }

    
    
    # Now count the number of each type
    foreach i $ilist {
	set v [lindex $plist [lsearch $itypelist $i] 1 ]
	if { [lsearch $itypelist $i] == -1 } {
	    mmsg::err [namespace current] "could not calculate proportions"
	}
	lset plist [lsearch $itypelist $i] 1 [expr $v + 1]
    }
    return $plist
}

# ::mbtools::utils::average --
#
# calculate the average of a list of numerical data
#
proc ::mbtools::utils::average { data { from 0 } { to 0 } } {
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

# ::mbtools::utils::stdev --
#
# calculate the stdev of a list of numerical data
#
proc ::mbtools::utils::stdev { data { from 0 } { to 0 } } {
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


# ::mbtools::utils::acorr --
#
# calculate an autocorrelation function from a list of numbers
#
proc ::mbtools::utils::acorr { data } {
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


# ::mbtools::utils::distance --
#
# Distance between two vectors
#
proc ::mbtools::utils::distance { vec1 vec2 } {


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


# ::mbtools::utils::normalize --
#
# Normalize a vector
#
proc ::mbtools::utils::normalize { vec } {
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

# ::mbtools::utils::scalevec --
#
# Scale a vector
#
proc ::mbtools::utils::scalevec { vec scale } {
    foreach val $vec {
	lappend svec [expr $val*$scale]
    }
    return $svec
}


#
#
#
# ::mbtools::utils::uniquelist --
#   
# Take a list of integers and construct a list that contains no
# duplication
#
#
proc ::mbtools::utils::uniquelist { original } { 


    lappend short [lindex $original 0]

    foreach element $original {
	if { [lsearch -exact $short $element] == -1 } {
	    lappend short $element
	}
    }

    return $short
}