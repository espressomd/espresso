puts "Program Information: \n[code_info]\n"
set starttime [expr int([lindex [exec date +%s])]];  

#############################################################
#  Convert to hexagonal lattice                             #
#############################################################
proc hexconvert { vec shift_vec d_space} {
    set hvec {{1. 0.5 0.5} {0.0 0.8660254 0.28867513} {0. 0. 0.81649658}}
    set rvec {0 0 0}
    set dim 3
    for {set j 0} { $j < $dim } {incr j} {
        for {set i 0} { $i < $dim } {incr i} {
             lset rvec $j [expr [lindex $rvec $j] + [lindex $vec $i] * [lindex $hvec $j $i]]
        }
        #lset rvec $j [expr ([lindex $d_space $j] * [lindex $rvec $j] + [lindex $shift_vec $j])]
    }
    lsqr $rvec
    for {set j 0} { $j < $dim } {incr j} {
        lset rvec $j [expr ([lindex $d_space $j] * [lindex $rvec $j] + [lindex $shift_vec $j])]
    }
    return $rvec
}

#############################################################
#  Create a single chain lattice                            #
#############################################################
#creates a linear chain starting from $part_id as id.
# $part_id: the first member of the chain
# $rvec: position of te first member
# $p_length: length of the chain
# $b_length: bond length between monomers
# this program returns the current number of particles 
proc create_chain { part_id rvec p_length b_length} {
    set posx [lindex $rvec 0]
    set posy [lindex $rvec 1]
    set posz [lindex $rvec 2]
    for {set j 0} { $j < $p_length } {incr j} {
        part $part_id pos $posx $posy $posz
        incr part_id
        set posz [expr $posz + $b_length]
    }
    return $part_id
}

# System identification: 
set name  "HMC"
set ident ""
set mypi          3.141592653589793
