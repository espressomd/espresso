# Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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

source "tests_common.tcl"

require_feature "INTER_DPD"
require_feature "TRANS_DPD"
require_feature "CONSTRAINTS"

puts "---------------------------------------------"
puts "- Testcase dpd.tcl running on [format %02d [setmd n_nodes]] nodes  -"
puts "---------------------------------------------"

if { [ catch {
# integration paramters
##################################################################
set n_part 11
set box_l 10
setmd box $box_l $box_l $box_l
setmd time_step 0.1
set dpd_gamma 10
set dpd_rcut 2
setmd periodic 1 1 1
setmd skin 0.0

# set up constraint
##################################################################
constraint wall normal 0 0 1 dist 0 type 0


# set up particles
##################################################################
for {set i 0} {$i < $n_part} {incr i} {
    set x 0
    set y 0
    set z [expr 0.2*($i+1)]

    part $i pos $x $y $z type 1 v 0 0 1
}


# set up dpd interaction
##################################################################
inter 0 1 inter_dpd $dpd_gamma $dpd_rcut 0 0 0 0

# calculate interactions
##################################################################
integrate 0


# verify results
##################################################################
for {set i 0} {$i < $n_part} {incr i} {
    #actual force
    set f_sim [part $i print f]

    #theoretical force
    set r "0 0 [lindex [part $i print pos] 2]"
    set v [part $i print v]
    set abs_r [expr sqrt(pow([lindex $r 0],2)+pow([lindex $r 1],2)+pow([lindex $r 2],2))]

    if {$abs_r > $dpd_rcut} {
	set f_theo "0 0 0"
    } else {
	set r_hat "[expr [lindex $r 0]/$abs_r] [expr [lindex $r 1]/$abs_r] [expr [lindex $r 2]/$abs_r]"
	set v_r "[expr [lindex $r_hat 0]*[lindex $v 0]+[lindex $r_hat 1]*[lindex $v 1]+[lindex $r_hat 2]*[lindex $v 2]]"
	set f_pre [expr -$dpd_gamma*pow((1-$abs_r/$dpd_rcut),2)*$v_r]
	set f_theo "[expr $f_pre*[lindex $r_hat 0]] [expr $f_pre*[lindex $r_hat 1]] [expr $f_pre*[lindex $r_hat 2]]"
    }
    
    set diff "[expr sqrt(pow([lindex $f_sim 0]-[lindex $f_theo 0],2) + pow([lindex $f_sim 1]-[lindex $f_theo 1],2) + pow([lindex $f_sim 2]-[lindex $f_theo 2],2))]"

    if {$diff > 1e-12} {
	error "INTER_DPD force derivation too large"
    }
}

# check TRANS_DPD
###################################################################
inter 0 1 inter_dpd 0 0 0 $dpd_gamma $dpd_rcut 0

for {set i 0} {$i < $n_part} {incr i} {
    part $i v 1 0 0
}

integrate 0


for {set i 0} {$i < $n_part} {incr i} {
    #actual force
    set f_sim [part $i print f]

    #theoretical force
    set r "0 0 [lindex [part $i print pos] 2]"
    set v [part $i print v]
    set abs_r [expr sqrt(pow([lindex $r 0],2)+pow([lindex $r 1],2)+pow([lindex $r 2],2))]

    if {$abs_r > $dpd_rcut} {
	set f_theo "0 0 0"
    } else {
	set r_hat "[expr [lindex $r 0]/$abs_r] [expr [lindex $r 1]/$abs_r] [expr [lindex $r 2]/$abs_r]"
	set I_minus_rxr "[expr 1-([lindex $r_hat 1]*[lindex $r_hat 2]-[lindex $r_hat 2]*[lindex $r_hat 1])] [expr 1-([lindex $r_hat 2]*[lindex $r_hat 0]-[lindex $r_hat 0]*[lindex $r_hat 2])] [expr 1-([lindex $r_hat 0]*[lindex $r_hat 1]-[lindex $r_hat 1]*[lindex $r_hat 0])]"
	set v_r "[expr [lindex $I_minus_rxr 0]*[lindex $v 0]] [expr [lindex $I_minus_rxr 1]*[lindex $v 1]] [expr [lindex $I_minus_rxr 2]*[lindex $v 2]]"
	set f_pre [expr -$dpd_gamma*pow((1-$abs_r/$dpd_rcut),2)]
	set f_theo "[expr $f_pre*[lindex $v_r 0]] [expr $f_pre*[lindex $v_r 1]] [expr $f_pre*[lindex $v_r 2]]"
    }
    
    set diff "[expr sqrt(pow([lindex $f_sim 0]-[lindex $f_theo 0],2) + pow([lindex $f_sim 1]-[lindex $f_theo 1],2) + pow([lindex $f_sim 2]-[lindex $f_theo 2],2))]"


    if {$diff > 1e-12} {
	error "TRANS_DPD force derivation too large"
    }
}   


# check with particle instead a constaint
####################################################################
# TRANS_DPD
####################################################################
constraint delete

part $n_part pos 0 0 0 type 0

integrate 0


for {set i 0} {$i < $n_part} {incr i} {
    #actual force
    set f_sim [part $i print f]

    #theoretical force
    set r "0 0 [lindex [part $i print pos] 2]"
    set v [part $i print v]
    set abs_r [expr sqrt(pow([lindex $r 0],2)+pow([lindex $r 1],2)+pow([lindex $r 2],2))]

    if {$abs_r > $dpd_rcut} {
	set f_theo "0 0 0"
    } else {
	set r_hat "[expr [lindex $r 0]/$abs_r] [expr [lindex $r 1]/$abs_r] [expr [lindex $r 2]/$abs_r]"
	set I_minus_rxr "[expr 1-([lindex $r_hat 1]*[lindex $r_hat 2]-[lindex $r_hat 2]*[lindex $r_hat 1])] [expr 1-([lindex $r_hat 2]*[lindex $r_hat 0]-[lindex $r_hat 0]*[lindex $r_hat 2])] [expr 1-([lindex $r_hat 0]*[lindex $r_hat 1]-[lindex $r_hat 1]*[lindex $r_hat 0])]"
	set v_r "[expr [lindex $I_minus_rxr 0]*[lindex $v 0]] [expr [lindex $I_minus_rxr 1]*[lindex $v 1]] [expr [lindex $I_minus_rxr 2]*[lindex $v 2]]"
	set f_pre [expr -$dpd_gamma*pow((1-$abs_r/$dpd_rcut),2)]
	set f_theo "[expr $f_pre*[lindex $v_r 0]] [expr $f_pre*[lindex $v_r 1]] [expr $f_pre*[lindex $v_r 2]]"
    }
    
    set diff "[expr sqrt(pow([lindex $f_sim 0]-[lindex $f_theo 0],2) + pow([lindex $f_sim 1]-[lindex $f_theo 1],2) + pow([lindex $f_sim 2]-[lindex $f_theo 2],2))]"



    if {$diff > 1e-12} {
	error "TRANS_DPD force derivation too large"
    }
}   

# INTER_DPD
###################################################################
inter 0 1 inter_dpd $dpd_gamma $dpd_rcut 0 0 0 0

for {set i 0} {$i < $n_part} {incr i} {
    part $i v 0 0 1
}


integrate 0

for {set i 0} {$i < $n_part} {incr i} {
    #actual force
    set f_sim [part $i print f]

    #theoretical force
    set r "0 0 [lindex [part $i print pos] 2]"
    set v [part $i print v]
    set abs_r [expr sqrt(pow([lindex $r 0],2)+pow([lindex $r 1],2)+pow([lindex $r 2],2))]

    if {$abs_r > $dpd_rcut} {
	set f_theo "0 0 0"
    } else {
	set r_hat "[expr [lindex $r 0]/$abs_r] [expr [lindex $r 1]/$abs_r] [expr [lindex $r 2]/$abs_r]"
	set v_r "[expr [lindex $r_hat 0]*[lindex $v 0]+[lindex $r_hat 1]*[lindex $v 1]+[lindex $r_hat 2]*[lindex $v 2]]"
	set f_pre [expr -$dpd_gamma*pow((1-$abs_r/$dpd_rcut),2)*$v_r]
	set f_theo "[expr $f_pre*[lindex $r_hat 0]] [expr $f_pre*[lindex $r_hat 1]] [expr $f_pre*[lindex $r_hat 2]]"
    }
    
    set diff "[expr sqrt(pow([lindex $f_sim 0]-[lindex $f_theo 0],2) + pow([lindex $f_sim 1]-[lindex $f_theo 1],2) + pow([lindex $f_sim 2]-[lindex $f_theo 2],2))]"


    if {$diff > 1e-12} {
	error "INTER_DPD force derivation too large"
    }
}

} res ] } {
    error_exit $res
}
exit 0
