#  This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
#  It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
#  and by which you are legally bound while utilizing this file in any form or way.
#  There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#  You should have received a copy of that license along with this program;
#  if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
#  write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
#  Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#  
#############################################################
#                                                           #
# ABHmath.tcl                                               #
# ===========                                               #
#                                                           #
# Script to supply some auxiliary mathematical functions.   #
#                                                           #
# Created:       24.02.2003 by BAM                          #
#                                                           #
#############################################################



proc PI { } {
# PI
# Returns an approximate value of the mathematical constant pi.
    return 3.141592653589793
}

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


proc pair_dist { part_id1 part_id2 } {
# pair_dist <part_id1> <part_id2> 
# Returns the distance of two particles with identities <part_id1> and <part_id2>.
    set pos1 [part $part_id1 print pos]    
    set pos2 [part $part_id2 print pos]
    set dist2 [expr sqr([lindex $pos1 0]- [lindex $pos2 0])]    
    set dist2 [expr $dist2 + sqr([lindex $pos1 1]- [lindex $pos2 1])]
    set dist2 [expr $dist2 + sqr([lindex $pos1 2]- [lindex $pos2 1])]
    return [expr sqrt($dist2)]
}

proc veclen {v} {
    set len [llength $v]
    set w [lindex $v 0]
    set res [expr $w*$w]
    for {set i 1} {$i < $len} {incr i} {
	set w [lindex $v $i]
	set res [expr $res + $w*$w]
    }
    return $res
}

proc vecadd {a b} {
    set len [llength $a]
    set res [expr [lindex $a 0] + [lindex $b 0]]
    for {set i 1} {$i < $len} {incr i} {
	set res "$res [expr [lindex $a $i] + [lindex $b $i]]"
    }
    return $res
}

proc vecsub {a b} {
    set len [llength $a]
    set res [expr [lindex $a 0] - [lindex $b 0]]
    for {set i 1} {$i < $len} {incr i} {
	set res "$res [expr [lindex $a $i] - [lindex $b $i]]"
    }
    return $res
}

proc vecscale {a v} {
    set len [llength $v]
    set res [expr $a*[lindex $v 0]]
    for {set i 1} {$i < $len} {incr i} {
	set res "$res [expr $a*[lindex $v $i]]"
    }
    return $res
}


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

    set a [expr $ssxy/$ssxx]
    set b [expr $avgy - $avgx*$a]

    set chi2 [expr ($ssyy - $a*$ssxy)/double($num-2)]
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
