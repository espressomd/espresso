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

proc LinRegression {l} {
    # l is a list {{x1 y1} {x2 y2} ...} of points.
    # the LinRegression returns the least-square linear fit a*x+b
    # and the standard errors da and db
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

    set tmp [expr ($ssyy - $a*$ssxy)/double($num-2)]
    if { $tmp < 0 } {puts "LinReg warning: s^2=$tmp should be > 0, set to 0"; set tmp 0}
    set s [expr sqrt($tmp)]
    set da [expr $s*sqrt(1.0/$num + $avgx*$avgx/$ssxx)]
    set db [expr $s/sqrt($ssxx)]

    return "$a $b $da $db"
}

