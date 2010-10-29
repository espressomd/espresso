#############################################################
#                                                           #
#  Template for Polyelectrloyte-Bundle Simulations          #
#                                                           #
#                                                           #
#  Created:       06.01.2004 by Sayar                       #
#  Last modified: 06.01.2004 by Sayar                       #
#                                                           #
#############################################################

puts " "
puts "======================================================="
puts "=                Polyelectrloyte-bundle_temp.tcl      ="
puts "======================================================="
puts " "

# System parameters
#############################################################
set f [open "input.txt" "r"]
while { [eof $f] == 0 } {
    gets $f temp;
    if {[lindex $temp 0] == "number_of_molecules" } { set n_poly [lindex $temp 1]} \
	elseif {[lindex $temp 0] == "polymer_length"} { set p_length [lindex $temp 1]} \
	elseif {[lindex $temp 0] == "c_ion_val"} { set c_ion_val [lindex $temp 1]} \
	elseif {[lindex $temp 0] == "density"} { set density [lindex $temp 1]} \
	elseif {[lindex $temp 0] == "cion_d"} { set cion_d [lindex $temp 1]} \
	elseif {[lindex $temp 0] == "int_steps"} { set int_steps [lindex $temp 1]} \
	elseif {[lindex $temp 0] == "int_n_times"} { set int_n_times [lindex $temp 1]} \
	elseif {[lindex $temp 0] == "ensemble_num"} { set ensemble_num [lindex $temp 1]} \
	elseif {[lindex $temp 0] == "fac_lb"} { set fac_lb [lindex $temp 1]} \
	elseif {[lindex $temp 0] == "loops"} { set loops [lindex $temp 1]} \
	elseif {[lindex $temp 0] == "runtime"} { set runtime [lindex $temp 1]} \
	elseif {[lindex $temp 0] == "lb_add"} { set lb_add [lindex $temp 1]} \
	elseif {[lindex $temp 0] == "bjerrum_list"} { set bjerrum_list [lrange $temp 1 end]}
}
close $f

set b_length 1.00
set d_max 2.24
set d_space [list [expr $d_max] [expr $d_max] 1]
#polymer backbone charge
set p_val 1 

# Interaction parameters
#############################################################

set ljr_cut       1.12246204831
set ljr_eps       1.0

set fene_r        2.0
set fene_k        7.0

set bend_k        100.0

set accuracy      1.0e-4

# Integration parameters
#############################################################

set time_step    0.005
set skin         0.5

set warm_steps   100
set warm_n_times 1

set min_dist     1.1

set ci_steps     100
set ci_n_times   1

set poly_steps     100
set poly_n_times   1

set ensequi_steps 100
set ensequi_n_times 1


# Other parameters
#############################################################
set tcl_precision 10
set mypi          3.141592653589793
set vmd_output    "no"
set vmd_wait      0

#############################################################
#  Setup System                                             #
#############################################################

set n_part [expr $n_poly * $p_length * (1. + 1./$c_ion_val) ]
set volume [expr $n_part/$density]
set sphere_rad [expr pow((3.0*$volume)/(4.0*$mypi),1.0/3.0)]
set  box_l       [expr 4.0*$sphere_rad + 6.0*$skin]
set box_l [expr pow($volume,1.0/3.0)]
set center      [expr $box_l/2.0]
set shift_vec [list [expr $box_l/2.0] [expr $box_l/2.0] [expr $box_l/2.0- ($p_length-1.)/2.*$b_length]]

setmd box_l     $box_l $box_l $box_l
setmd periodic  1 1 1 
setmd time_step $time_step
setmd skin      $skin
setmd gamma     1.0
setmd temp      1.0

# Interaction setup
#############################################################

# type 0 charged backbone,  type 1 counterion
# repulsive LJ for all
set ljr_shift  0.25
inter 0 0 lennard-jones $ljr_eps 1.0 $ljr_cut $ljr_shift 0
inter 0 1 lennard-jones $ljr_eps $cion_d $ljr_cut $ljr_shift 0
inter 1 1 lennard-jones $ljr_eps $cion_d $ljr_cut $ljr_shift 0

# FENE 
inter 0 fene $fene_k $fene_r

# Stiffness
inter 1 angle $bend_k

#############################################################
#  Create a bundle of n_poly chains                         #
#############################################################

set part_id 0
set input_file [open "hex_input.txt" "r"]
for {set i 0} { $i < [expr $n_poly] } {incr i} {
    gets $input_file vec
    set rvec [hexconvert $vec $shift_vec $d_space]
    set part_id [create_chain $part_id $rvec $p_length $b_length]
}
close $input_file
set input_file [open "hex_input.txt" "r"]
set shift_vec [list [expr $box_l/2.0+1.] [expr $box_l/2.0+1.] [expr $box_l/2.0- ($p_length-1.)/2.*$b_length]]
for {set i 0} { $i < [expr $n_poly] } {incr i} {
    gets $input_file vec
    set rvec [hexconvert $vec $shift_vec $d_space]
    set part_id [create_chain $part_id $rvec [expr $p_length/$c_ion_val] [expr $b_length*$c_ion_val]]
}
close $input_file

#fix the particles on the chains and set up lj types and add stiffness and put the bonds
for {set i 0} { $i < [expr $n_poly * $p_length ]} {incr i} {
    part $i type 0 q $p_val 
    if { [expr $i % $p_length ]  != 0 } {
        part [expr $i-1] bond 0 $i
        if { [expr $i % $p_length]  != 1 } {
            part [expr $i-1] bond 1 [expr $i-2] [expr $i]
        }
    }
}

# particle numbers
set n_mono [expr $n_poly * $p_length]
set n_ci   [expr $n_mono/$c_ion_val]

# Counterions 
for {set i 0} { $i < $n_ci } {incr i} {
    part [expr $n_mono + $i]  type 1 q -$c_ion_val
}

# Status report
puts "Simulate $n_part particles in a spherical cell with radius $sphere_rad"
puts "$n_poly PE backbones of length $p_length and charge distance $b_length"
puts "Each monomer has $p_val charge"
puts "neutralized by $n_ci counterions with valency $c_ion_val."
puts "Interactions:\n[inter]"

#############################################################
#  VMD connection                                           #
#############################################################
if { $vmd_output=="yes" } {
    prepare_vmd_connection "$name$ident" $vmd_wait 1
}
