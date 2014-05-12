# Copyright (C) 2010,2012,2013 The ESPResSo project
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
    namespace export min_vec 
    namespace export distance
    namespace export distance_min
    namespace export find_proportions
    namespace export matrix_multiply
    namespace export matrix_vec_multiply
    namespace export cross_product
    namespace export rotation_matrix
    namespace export add_vecs
    namespace export dot_product
    namespace export perp_vec
}


# ::mbtools::utils::dot_product
proc ::mbtools::utils::dot_product { A B } {
    set dp 0
    for { set i 0 } { $i < 3 } { incr i } {
	set dp [expr $dp + [lindex $A $i]*[lindex $B $i] ]
    }
    return $dp

}

# ::mbtools::utils::perp_vec
#
#  returns a vector perpendicular to A
#  nothing special about the vector except that it's one of the perpendicular options and is normalized
#
proc ::mbtools::utils::perp_vec {A} {
    set x [lindex $A 0]
    set y [lindex $A 1]
    if {$x  == 0} {
	set x 1
    } elseif {$y == 0} {
	set y 1
    } else {
	set ratio [expr $y / $x]
	set y [expr $x * $ratio * 2]
    }
    set some_vector "$x $y [lindex $A 2]"
    set some_vector [::mbtools::utils::normalize $some_vector]
    return [::mbtools::utils::cross_product $A $some_vector]
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

# ::mbtools::utils::min_vec --
#
# Minimal vector between two vectors
#
proc ::mbtools::utils::min_vec { vec1 vec2 } {


    if { [ llength $vec1 ] != [llength $vec2] } {
	mmsg::err [namespace current] "cannot find distance between vectors $vec1 and $vec2 because they are different lengths"
    }
    
    set box [setmd box_l]
    set sum 0
    for {set i 0 } { $i < [llength $vec1] } { incr i } {
	lappend diff [expr [lindex $vec1 $i] - [lindex $vec2 $i]]
	lappend add [expr -floor(([lindex $diff $i]/[lindex $box $i])+0.5)*[lindex $box $i]]
    }
    set vec [vecadd $diff $add]
    return $vec
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


# ::mbtools::utils::distance_min --
#
# Minimal distance between two vectors
#
proc ::mbtools::utils::distance_min { vec1 vec2 } {


    if { [ llength $vec1 ] != [llength $vec2] } {
	mmsg::err [namespace current] "cannot find distance between vectors $vec1 and $vec2 because they are different lengths"
    }
    
    set box [setmd box_l]
    set sum 0
    for {set i 0 } { $i < [llength $vec1] } { incr i } {
	lappend diff [expr [lindex $vec1 $i] - [lindex $vec2 $i]]
	lappend add [expr -floor(([lindex $diff $i]/[lindex $box $i])+0.5)*[lindex $box $i]]
    }
    set vec [vecadd $diff $add]
    return [veclen $vec]
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

#
#
# :: mbtools::utils::rotate --
#
# Take a {x y z} point and rotate it an angle sigma around either
# the x y or z axis.
#
#
proc ::mbtools::utils::rotation_matrix {axis phi} {
    
    set x [lindex $axis 0]
    set y [lindex $axis 1]
    set z [lindex $axis 2]
    set x2 [expr $x * $x]
    set y2 [expr $y * $y]
    set z2 [expr $z * $z]

    set rcos  [expr cos($phi)]
    set rsin  [expr sin($phi)]

    set row0 ""
    set row1 ""
    set row2 ""

    lappend row0 [expr $x2 + ($y2+$z2)*$rcos]
    lappend row0 [expr $x*$y*(1.0-$rcos) - $z*$rsin]
    lappend row0 [expr $x*$z*(1.0-$rcos) + $y*$rsin]
    lappend row1 [expr $x*$y*(1.0-$rcos) + $z*$rsin]
    lappend row1 [expr $y2 + ($x2+$z2)*$rcos]
    lappend row1 [expr $y*$z*(1.0-$rcos) - $x*$rsin]
    lappend row2 [expr $x*$z*(1.0-$rcos) - $y*$rsin]
    lappend row2 [expr $y*$z*(1.0-$rcos) + $x*$rsin]
    lappend row2 [expr $z2 + ($x2+$y2)*$rcos]


    set rotation_matrix ""
    lappend rotation_matrix $row0
    lappend rotation_matrix $row1
    lappend rotation_matrix $row2

    return $rotation_matrix
}

# ::mbtools::utils::cross_product
proc ::mbtools::utils::cross_product { A B } {
    set cp ""
    set Ax [lindex $A 0]
    set Ay [lindex $A 1]
    set Az [lindex $A 2]
    set Bx [lindex $B 0]
    set By [lindex $B 1]
    set Bz [lindex $B 2]

    lappend cp [expr  $Ay*$Bz - $Az*$By]
    lappend cp [expr -$Ax*$Bz + $Az*$Bx]
    lappend cp [expr  $Ax*$By - $Ay*$Bx]

    return $cp

}

# ::mbtools::utils::add_vecs
proc ::mbtools::utils::add_vecs { A B} {
    lappend C [expr [lindex $A 0] + [lindex $B 0]]    
    lappend C [expr [lindex $A 1] + [lindex $B 1]]
    lappend C [expr [lindex $A 2] + [lindex $B 2]]
    return $C
}

# ::mbtools::utils::matrix_multiply
proc ::mbtools::utils::matrix_multiply {A B} {
    set C ""
    for {set i 0} {$i < 3} {incr i} {
	set row ""
	for {set j 0} {$j < 3} {incr j} {
	set sum 0
	    for {set k 0} {$k < 3} {incr k} {
		set sum [expr $sum + [lindex [lindex $A $i] $k] * [lindex [lindex $B $k] $j]]
	    }
	    lappend row $sum
	}
	lappend C $row
    }
    return $C
}
