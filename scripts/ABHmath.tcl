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
