#############################################################
#                                                           #
# ABHmath.tcl                                               #
# ===========                                               #
#                                                           #
# Script to supply some auxiliary mathematical functions.   #
#                                                           #
#############################################################
#
# Copyright (C) 2010,2012,2013 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#   Max-Planck-Institute for Polymer Research, Theory Group
#  
# This file is part of ESPResSo.
#  
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 
#  

#############################################################
# CONSTANTS
#############################################################

proc PI { } {
# PI
# Returns an approximate value of the mathematical constant pi.
    return 3.141592653589793
}

proc KBOLTZ { } {
    # Returns Boltzmann constant in Joule/Kelvin
    return 1.380658e-23
}

proc ECHARGE { } {
    # Returns elementary charge in Coulomb
    return 1.60217733e-19
}

proc NAVOGADRO { } {
    # Returns Avogadro number
    return 6.0221367e23
}

proc SPEEDOFLIGHT { } {
    # Returns speed of light in meter/second
    return 2.99792458e8
}

proc EPSILON0 { } {
    # Returns dielectric constant of vaccum in Coulomb^2/(Joule meter)
    return 8.854187817e-12
}

proc ATOMICMASS { } {
    # Returns the atomic mass unit u in kilogramms
    return 1.66053873e-27
}

#############################################################
# MATHEMATICAL FUNCTIONS
#############################################################

proc sqr { arg } {
# sqr <arg>
# Returns the square of <arg>.
    return [expr $arg*$arg]
}

proc min { arg1 arg2 } {
# min <arg1> <arg2>
# Returns the minimum of <arg1> and <arg2>.
    if { $arg1 < $arg2 } { return $arg1 } else { return $arg2 }
}

proc max { arg1 arg2 } {
# max <arg1> <arg2>
# Returns the maximum of <arg1> and <arg2>.
    if { $arg1 > $arg2 } { return $arg1 } else { return $arg2 }
}

proc sign { arg } {
# returns the signum-function of <arg>
    if { $arg==0 } { return $arg } { return [expr round($arg/abs($arg))] }
}

#############################################################
# RANDOM FUNCTIONS
#############################################################

proc gauss_random { } {
    # returns random numbers which have a Gaussian distribution
    set v1  [expr 2.0*[t_random]-1.0]
    set v2  [expr 2.0*[t_random]-1.0]
    set rsq [expr [sqr $v1]+[sqr $v2]]
    while { $rsq >= 1.0 || $rsq == 0.0 } {
	set v1  [expr 2.0*[t_random]-1.0]
	set v2  [expr 2.0*[t_random]-1.0]
	set rsq [expr [sqr $v1]+[sqr $v2]]
    }
    set fac [expr sqrt(-2.0*log($rsq)/$rsq)]
    return [expr $v2*$fac]
}

proc dist_random { dist { max "1.0" } } {
    # returns random numbers in the interval [0,1] which have a distribution 
    # according to the distribution function p(x) <dist> which has to be given 
    # as a tcl list containing equally spaced values of p(x). 
    # If p(x) contains values larger than 1 the maximum or any number larger 
    # than that has to be given <max>

    set bins [expr [llength $dist] -1]

    while { 1 } {
	set v1 [t_random]
	set v2 [expr $max*[t_random]]
	set ind [expr int($v1*$bins)]
	set arg [expr ($v1*$bins)-$ind]
	if { $v2 < [expr [lindex $dist $ind] + ([lindex $dist [expr $ind+1]]-[lindex $dist $ind])*$arg]} { break }
    }
    return $v1
}

proc vec_random { {len "1.0"} } {
    # returns a random vector of length len
    # (uniform distribution on a sphere)
    # This is done by chosing 3 uniformly distributed random numbers [-1,1]
    # If the length of the resulting vector is <= 1.0 the vector is taken and normalized
    # to the desired length, otherwise the procedure is repeated until succes.
    # On average the procedure needs 5.739 random numbers per vector.
    # (This is probably not the most efficient way, but it works!)
    # Ask your favorit mathematician for a proof!
    while {1} {
        for {set i 0} {$i<3} {incr i} { lappend vec [expr 2.0*[t_random]-1.0] }
        if { [veclen $vec] <= 1.0 } { break }
        unset vec
    }
    return [vecnorm $vec $len]
}

proc phivec_random { v phi {len "1.0"} } {
    # return a random vector at angle phi with v and length len
    # find orthonormal basis
    set v  [vecnorm $v] 
    set bas1 [orthovec3d $v]
    set bas2 [veccross_product3d $v $bas1]
    # select random angle theta
    set theta [expr [t_random]*[PI]*2.0]
    # build vector
    set res [vecnorm $v [expr -cos($phi)]]
    set res [vecadd $res [vecnorm $bas1 [expr sin($phi)*cos($theta)]] ]
    set res [vecadd $res [vecnorm $bas2 [expr sin($phi)*sin($theta)]] ]
    # normalize to length
    return [vecnorm $res $len]
}


#############################################################
# PARTICLE OPERATIONS
#############################################################

# Operations involving particle positions.
# The parameters pi can either denote the particle identity 
# (then the particle position is extracted with the part command) 
# or the particle position directly
# When the optional box parameter for minimum image conventions is 
# omited the functions use the setmd box_l command.

proc bond_vec { p1 p2 } {
    # Calculate bond vector pointing from  particles p2 to p1
    # return = (p1.pos - p2.pos)
    if { [llength $p1]==1 } { set p1 [part $p1 print pos] }
    if { [llength $p2]==1 } { set p2 [part $p2 print pos] }
    return [vecsub $p1 $p2]
}

proc bond_vec_min { p1 p2 {box ""} } {
    # Calculate bond vector pointing from  particles p2 to p1
    # return = (p1.pos - p2.pos)
    if { [llength $p1]==1 } { set p1 [part $p1 print pos] }
    if { [llength $p2]==1 } { set p2 [part $p2 print pos] }
    if { $box == "" }       { set box [setmd box_l] }
    set vec [vecsub $p1 $p2]
    # apply minimal image convention
    foreach v $vec b $box { lappend add [expr -floor(($v/$b)+0.5)*$b] }
    return [vecadd $vec $add]
}

proc bond_length { p1 p2 } {
    # Calculate bond length between particles p1 and p2
    if { [llength $p1]==1 } { set p1 [part $p1 print pos] }
    if { [llength $p2]==1 } { set p2 [part $p2 print pos] }
    set vec [vecsub $p1 $p2]
    return [veclen $vec]
}

proc bond_length_min { p1 p2 {box ""} } {
    # Calculate minimum image bond length between particles p1 and p2
    if { [llength $p1]==1 } { set p1 [part $p1 print pos] }
    if { [llength $p2]==1 } { set p2 [part $p2 print pos] }
    if { $box == "" }       { set box [setmd box_l] }
    set vec [vecsub $p1 $p2]
    # apply minimal image convention
    foreach v $vec b $box { lappend add [expr -floor(($v/$b)+0.5)*$b] }
    set vec [vecadd $vec $add]
    return [veclen $vec]
}

proc bond_angle { p1 p2 p3 {type "r"} } {
    # Calculate bond angle between particles p1, p2 and p3
    if { [llength $p1]==1 } { set p1 [part $p1 print pos] }
    if { [llength $p2]==1 } { set p2 [part $p2 print pos] }
    if { [llength $p3]==1 } { set p3 [part $p3 print pos] }
    set vec1 [unitvec $p1 $p2]
    set vec2 [unitvec $p3 $p2]
    set cosphi [vecdot_product $vec1 $vec2]
    if { $cosphi<=-1 || $cosphi >= 1} { set phi [PI]  } else {
	set phi  [expr acos($cosphi)] }
    if { $type == "d" } { return [expr 180.0*$phi/[PI]] } 
    return $phi 
}

proc bond_angle_min { p1 p2 p3 {type "r"} } {
    # Calculate bond angle between particles p1, p2 and p3 using the minimum image convention
    if { [llength $p1]==1 } { set p1 [part $p1 print pos] }
    if { [llength $p2]==1 } { set p2 [part $p2 print pos] }
    if { [llength $p3]==1 } { set p3 [part $p3 print pos] }
    set vec1 [bond_vec_min $p1 $p2]
    set vec2 [bond_vec_min $p3 $p2]
    set phi  [expr acos([vecdot_product $vec1 $vec2] / ([veclen $vec1] * [veclen $vec2]))]
    if { $type == "d" } { return [expr 180.0*$phi/[PI]] } 
    return $phi 
}

proc bond_dihedral { p1 p2 p3 p4 {type "r"} } {
    # Calculate bond dihedral between particles p1, p2, p3 and p4
    # Type defines the angle convention and the unit
    # First letter: r - radiant; d - degree
    # second letter: "" - 0,360; b - bio convention -180,180
    if { [llength $p1]==1 } { set p1 [part $p1 print pos] }
    if { [llength $p2]==1 } { set p2 [part $p2 print pos] }
    if { [llength $p3]==1 } { set p3 [part $p3 print pos] }
    if { [llength $p4]==1 } { set p4 [part $p4 print pos] }

    set r_12 [unitvec $p2 $p1]
    set r_32 [unitvec $p2 $p3]
    set r_34 [unitvec $p4 $p3]

    set aXb  [vecnorm [veccross_product3d $r_12 $r_32]]
    set bXc  [vecnorm [veccross_product3d $r_32 $r_34]]

    set cosphi [vecdot_product $aXb $bXc]
    if { $cosphi<-1 || $cosphi > 1} { set phi [PI]  } else {
	set phi [expr acos($cosphi)] }

    if { [vecdot_product $aXb $r_34] < 0 } { set phi [expr 2*[PI]-$phi] }
    if { $type == "d" } { return [expr 180.0*$phi/[PI]] } 
    if { $type == "db" } { 
	if { $phi<[PI] } { return [expr 180.0*$phi/[PI]]
	} else {        return [expr (180.0*$phi/[PI])-360] } } 
    return $phi
}

proc part_at_dist { p dist } {
    # return position of a new particle at distance dist from p
    # with random orientation
    if { [llength $p]==1 } { set p [part $p print pos] }
    set vec [vec_random $dist]
    return [vecadd $p $vec]
}

proc part_at_angle { p1 p2 phi {len "1.0"} } {
    # return particle at distance len (default=1.0) from p2 which builds
    # a bond angle phi for (p1, p2, res) 
    if { [llength $p1]==1 } { set p1 [part $p1 print pos] }
    if { [llength $p2]==1 } { set p2 [part $p2 print pos] }
    set vec [ phivec_random [bond_vec $p2 $p1] $phi $len ]
    return [vecadd $p2 $vec]
}

proc part_at_dihedral { p1 p2 p3 theta {phi "rnd"} {len "1.0"} } {
    # return particle at distance len (default=1.0) from p3 which builds
    # a bond angle phi (default=random) for (p2, p3, res) and a
    # dihedral angle theta for (p1, p2, p3, res) 
    if { [llength $p1]==1 } { set p1 [part $p1 print pos] }
    if { [llength $p2]==1 } { set p2 [part $p2 print pos] }
    if { [llength $p3]==1 } { set p3 [part $p3 print pos] }
    set v12 [bond_vec $p2 $p1]
    set v23 [bond_vec $p3 $p2]
    set vec [create_dihedral_vec $v12 $v23 $theta $phi $len]
    return [vecadd $p3 $vec]
}

#############################################################
# INTERACTION RELATED
#############################################################

# Calculate the shift value for the lennard jones interaction
# This function is obsolete. Better use a value of "auto" for the
# shift using "inter" to set up the potential.
proc calc_lj_shift { lj_sigma lj_cutoff } {
    set fcut [expr $lj_sigma / (1.0*$lj_cutoff) ]
    return [expr -(pow($fcut,12)-pow($fcut,6))]
}

#############################################################
# VECTOR OPERATIONS
#############################################################

# A vector is a tcl list of numbers
# Some functions are provided only for three dimensional vectors.
# corresponding functions contain 3d at the end of the name

proc veclen {v} {
    # return the length of a vector
    set res 0
    foreach e $v { set res [expr $res + ($e*$e)] } 
    return [expr sqrt($res)]
}

proc veclensqr {v} {
    # return the length of a vector squared
    set res 0
    foreach e $v { set res [expr $res + ($e*$e)] } 
    return $res
}

proc vecadd {a b} {
    # add vector a to vector b: return = (a+b)
    foreach ai $a bi $b { lappend res [expr $ai + $bi] }
    return $res
}

proc vecsub {a b} {
    # subtract vector b from vector a: return = (a-b)
    foreach ai $a bi $b { lappend res [expr $ai - $bi] }
    return $res
}

proc vecscale {s v} {
    # scale vector v with factor s: return = (s*v)
    foreach e $v { lappend res [expr $s * $e] } 
    return $res
}

proc vecdot_product { a b } {
    # calculate dot product of vectors a and b: return = (a.b)
    set res 0
    foreach ai $a bi $b { set res [expr $res + ($ai*$bi)] }
    return $res
}

proc veccross_product3d { a b } {
    # calculate the cross product of vectors a and b: return = (a x b)
    for {set i 0} {$i<3} {incr i} { 
	lappend res [expr [lindex $a [expr ($i+1)%3]]*[lindex $b [expr ($i+2)%3]] - [lindex $a [expr ($i+2)%3]]*[lindex $b [expr ($i+1)%3]] ]
    }
    return $res
}

proc matvec_product { m v } {
    # Caclulate the dot product m.v, where m is a (n x n) matrix 
    # and v is a vector with n elements
    set n [llength $m]; set res ""
    for { set i 0 } { $i < $n } { incr i } {
	set tmp 0
	for { set j 0 } { $j < $n } { incr j } { 
	    set tmp [expr $tmp+([lindex $m $i $j]*[lindex $v $j])]
	}
	lappend res $tmp
    }
    return $res
}

proc vecnorm { v {len "1.0"} } {
    # normalize a vector to length len (default 1.0) 
    set scale [expr $len/[veclen $v]]
    return [vecscale $scale $v]
}

proc unitvec { p1 p2 } {
    # return unit vector pointing from p1 to p2
    set res [vecsub $p2 $p1]
    return [vecnorm $res]
}

proc orthovec3d { v {len "1.0"} } {
    # return orthogonal vector to v with length len (default 1.0)
    # This vector does not have a random orientation in the plane perpendicular to v
    set nzero 0
    for { set i 0 } { $i < 3 } { incr i } {
	if { [lindex $v $i] == 0  } { 
	    incr nzero; set zind $i
	} else { set nzind $i }   
    }
    if { $nzero == 0 } {
	set res [lindex $v 0]
	lappend res [lindex $v 1]
	set tmp [expr [sqr [lindex $v 0]]+[sqr [lindex $v 1]]]
	lappend res [expr -$tmp/[lindex $v 2] ]
    }
    if { $nzero == 1 } { set res { 0 0 0 }; lset res $zind  1 }
    if { $nzero == 2 } { set res { 1 1 1 }; lset res $nzind 0 }
    if { $nzero == 3 } {
	puts "Can not handle null vector"
	return TCL_ERROR
    }
    return [vecnorm $res $len]
}

proc create_dihedral_vec { vec1 vec2 theta {phi "rnd"} {len "1.0"} } {
    # create last vector of a dihedral (vec1, vec2, res) with dihedral angle
    # theta and bond angle (vec2, res) phi and length len.
    # construct orthonormal basis (vec2, b, c)
    set vec1 [vecnorm $vec1]
    set vec2 [vecnorm $vec2]
    set a [vecscale [vecdot_product $vec1 $vec2] $vec2]
    set b [vecsub $a $vec1]
    set b [vecnorm $b]
    set c [veccross_product3d $vec2 $b]
    # place dihedral angle theta in plane
    set d [vecadd [vecscale [expr cos($theta)] $b] [vecscale [expr sin($theta)] $c] ]
    # choose bond angle phi
    if { $phi == "rnd" } { set phi [expr 1.0*[PI]*[t_random]] }
    set res [vecadd [vecscale [expr sin($phi)] $d] [vecscale [expr -cos($phi)] $vec2]]
    
    return [vecnorm $res $len]
}

#############################################################
# TCL LIST OPERATIONS
#############################################################

proc average {list_to_average} {
    # Returns the avarage of the provided List
    set avg 0.0
    foreach avg_i $list_to_average { set avg [expr $avg + $avg_i] }
    return [expr $avg/(1.0*[llength $list_to_average])]
}

proc list_add_value { list val } {
    # Add <val> to each element of <list>
    for { set i 0 } { $i<[llength $list] } { incr i } {
	lappend res [expr [lindex $list $i] + $val]
    }
    return $res
}

proc flatten { list } {
    # flattens a nested list
    set res ""
    foreach sl $list {
	if { [llength $sl] > 1 } { set sl [flatten $sl] } 
	foreach e $sl { lappend res $e }
    }
    return $res
}

# Checks wether list contains val. returns the number of occurences of val in list.
proc list_contains { list val } {
    set res 0
    foreach e $list { if { $e == $val } { incr res } }
    return $res
}


# nlreplace returns a new list formed by replacing the element 
# indicated by <ind> with the element <arguments>. <ind> is a tcl list of 
# numbers specifying an element in a nested list which must exist
# (similar to lindex).
# list - tcl list (may be nested)
# ind  - tcl list specifying an element in a nested list
# elemnt - R 
proc nlreplace { list ind element } {
    if { [llength $ind] > 1 } {
	set i   [lindex $ind 0]
	set tmp [lindex $list $i]
	set ind [lreplace $ind 0 0]
	set tmp [nlreplace $tmp $ind $element]
	set list [lreplace $list $i $i $tmp]
    } else {
	set list [lreplace $list $ind $ind $element]
    }
    return $list
}

#############################################################
# Histograms
#############################################################

proc histogram { list nbins { min "auto" } { max "auto" }  { type "lin" } { ret "xy" } } {
    # calculate a histogram for the values given in list
    # The histogram with nbins bins accounts for values within the range min to max
    # The bins are either aequidistant or logarythmically aequidistant
    # The histogram is returned in xy pairs or only a list of histogram values
    # There is no normalisation since it depends too much on what you want 
    set num [llength $list]
    if { $min == "auto" || $max == "auto" } {
	set min [lindex $list 0]; set max $min
	for { set i 1 } { $i < $num } { incr i } {
	    if { $min > [lindex $list $i] } { set min [lindex $list $i] }
	    if { $max < [lindex $list $i] } { set max [lindex $list $i] }
	} 
    } 

    if { $type == "lin" || $type == "linear" } {	
	# linear case
	set fac [expr $nbins/(1.0*$max-$min)]
	for  { set i 0 } { $i < $nbins } { incr i } { lappend res 0 
	    lappend xval [expr $min+(($i+0.5)/$fac)] }
	for { set i 1 } { $i < $num } { incr i } {
	    set ind [expr int(([lindex $list $i]-$min)*$fac)]
	    if { $ind >= 0 && $ind < $nbins } {  lset res $ind [expr [lindex $res $ind]+1] }
	}
   } else {
       # logarithmic case 
       set fac [expr $nbins/(1.0*log($max)-log($min))]
       set fac2 [expr pow(($max/$min),1.0/(1.0*$nbins))]
       set tmp [expr $min/(sqrt($fac2))]
       for  { set i 0 } { $i < $nbins } { incr i } { lappend res 0 
	   set tmp [expr $tmp*$fac2]
	   lappend xval $tmp }
       for { set i 0 } { $i < $num } { incr i } {
	   set ind [expr int((log([lindex $list $i])-log($min))*$fac)]
	   if { $ind >= 0 && $ind < $nbins } {  lset res $ind [expr [lindex $res $ind]+1] }
      }
    } 
    if { $ret == "xy" } { 
	for  { set i 0 } { $i < $nbins } { incr i } { lappend res2 "[lindex $xval $i] [lindex $res $i]" } 
	return $res2
    } else { return $res }
}

#############################################################
# REGRESSION
#############################################################

proc LinRegression {l} {
    # l is a list {{x1 y1} {x2 y2} ...} of points.
    # LinRegression returns the least-square linear fit a*x+b
    # and the standard errors da and db.

    set num [llength $l]

    if { $num <= 2 } {
	error "not enough points for regression"
    }

    set sumx 0.0
    set sumy 0.0
    set sumx2 0.0
    set sumy2 0.0
    set sumxy 0.0

    foreach data $l {
	set x [lindex $data 0]
        set y [lindex $data 1]

        set sumx [expr $x + $sumx]
        set sumx2 [expr $x*$x + $sumx2]
        set sumy [expr $y + $sumy]
        set sumy2 [expr $y*$y + $sumy2]
        set sumxy [expr $x*$y + $sumxy]
    }

    set avgx [expr $sumx/$num]
    set avgy [expr $sumy/$num]
    set ssxx [expr $sumx2 - $num*$avgx*$avgx]
    if { $ssxx < 1e-4 } {
	error "data points too close together"
    }
    set ssyy [expr $sumy2 - $num*$avgy*$avgy]
    set ssxy [expr $sumxy - $num*$avgx*$avgy]

    set b [expr $ssxy/$ssxx]
    set a [expr $avgy - $avgx*$b]

    set chi2 [expr ($ssyy - $b*$ssxy)/double($num-2)]
    if { $chi2 < 0 } {puts "LinReg warning: s^2=$chi2 should be > 0, set to 0"; set tmp 0}
    set chi [expr sqrt($chi2)]
    set da [expr $chi*sqrt(1.0/$num + $avgx*$avgx/$ssxx)]
    set db [expr $chi/sqrt($ssxx)]

    return "$a $b $da $db"
}

proc LinRegressionWithSigma {l} {
    # l is a list {{x1 y1 s1} {x2 y2 s2} ...} of points with standard deviations.
    # LinRegression returns the least-square linear fit a*x+b
    # plus the standard errors da and db, cov(a,b) and chi.
    set num [llength $l]

    if { $num <= 2 } {
	error "not enough points for regression"
    }

    set s 0.0
    set sx 0.0
    set sy 0.0
    set sxx 0.0
    set sxy 0.0

    foreach data $l {
	set x [lindex $data 0]
        set y [lindex $data 1]
        set sigma [lindex $data 2]

	set wt [expr 1/($sigma*$sigma)]
	set s   [expr $s   + $wt]
        set sx  [expr $sx  + $wt*$x]
        set sxx [expr $sxx + $wt*$x*$x]
        set sy  [expr $sy  + $wt*$y]
        set sxy [expr $sxy + $wt*$x*$y]
    }
    if { $s < 1e-4 } {
	error "weights of data points too small"
    }
    
    set si [expr 1/$s]
    set stt [expr $sxx - $sx*$sx*$si]

    set a [expr ($sxy - $sx*$sy*$si)/$stt]
    set b [expr ($sy - $a*$sx)*$si]

    set da [expr sqrt((1 + $sx*$sx/($s*$stt))*$si)]
    set db [expr sqrt(1/$stt)]

    set chi2 0
    foreach data $l {
	set x [lindex $data 0]
        set y [lindex $data 1]
        set sigma [lindex $data 2]
	set diff [expr ($y - $b - $a*$x)/$sigma]
	set chi2 [expr $chi2 + $diff*$diff]
    }
    set chi [expr sqrt($chi2)]
    set covab [expr -$sx/($s*$stt)]

    return "$a $b $da $db $covab $chi"
}




#############################################################
# Compatibility Issues
#############################################################

# !!!!!!!!!!
# DO NOT USE THE FOLLOWING PROCEDURES ANYMORE
# THEY ARE ONLY HERE DUE TO BACKWARDS COMPATIBILITY
# !!!!!!!!!!

proc lsqr { v } {
    # return the length of a vector squared
    set res 0
    foreach e $v { set res [expr $res + ($e*$e)] } 
    return $res
}

proc dot_product { a b } {
    # calculate dot product of vectors a and b: return = (a.b)
    set res 0
    foreach ai $a bi $b { set res [expr $res + ($ai*$bi)] }
    return $res
}

proc find_unit_vector { vec1 vec2 } {
    return [unitvec $vec1 $vec2]
}

proc pair_vec { part_id1 part_id2 } {
    return [bond_vec $part_id1 $part_id2]
}

proc pair_dist { part_id1 part_id2 } {
    # pair_dist <part_id1> <part_id2> 
    # Returns the distance between two particles with identities <part_id1> and <part_id2>.
    if { [llength $part_id1]==1 } { set part_id1 [part $part_id1 print pos] }
    if { [llength $part_id2]==1 } { set part_id2 [part $part_id2 print pos] }
    set vec [vecsub $part_id1 $part_id2]
    return [veclen $vec]
}

proc find_vector_mag { vec1 vec2} {
    # Definition of these two functions is opposite
    return [bond_vec $vec2 $vec1]
}

proc g_random { } {
    # returns random numbers which have a Gaussian distribution
    set v1  [expr 2.0*[t_random]-1.0]
    set v2  [expr 2.0*[t_random]-1.0]
    set rsq [expr [sqr $v1]+[sqr $v2]]
    while { $rsq >= 1.0 || $rsq == 0.0 } {
	set v1  [expr 2.0*[t_random]-1.0]
	set v2  [expr 2.0*[t_random]-1.0]
	set rsq [expr [sqr $v1]+[sqr $v2]]
    }
    set fac [expr sqrt(-2.0*log($rsq)/$rsq)]
    return [expr $v2*$fac]
}

#############################################################
#
# calculate minimal image length of the vector connecting 
# 2 points inside periodic box box (slow and dirty!)
#
#############################################################
proc min_img_bond_length { pos1 pos2 box} {
    set dim [llength $pos1]
    set res 0.0
    for {set j 0} { $j < $dim } {incr j} {
	set dist [expr abs([lindex $pos1 $j]-[lindex $pos2 $j])]
	while { $dist > [expr [lindex $box $j]/2.0] } { 
	    set dist [expr $dist - [lindex $box $j]] 
	}
	set res [expr $res+ [sqr $dist] ]
    }
    return [expr sqrt($res)]
}
