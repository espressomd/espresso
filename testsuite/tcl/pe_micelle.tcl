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

#############################################################
#  calculate square of a vector                             #
#############################################################
proc lsqr { arg } {
# sqr <arg>
# Returns the square of <arg>.
    set dim [llength $arg]
    set sum 0.0
    for {set xi 0} {$xi < $dim} {incr xi} {
        set sum [expr $sum + [lindex $arg $xi] * [lindex $arg $xi] ]
    }
    #puts "the vector $arg has length^2 $sum"
    return $sum
}


#############################################################
#  Convert to hexagonal lattice                             #
#############################################################
proc hexconvert { vec shift_vec d_space} {
    set hvec {{1. 0.5 0.5} {0.0 0.8660254 0.28867513} {0. 0. 0.81649658}}
    set rvec {0 0 0}
    set dim 3
    for {set j 0} { $j < $dim } {incr j} {
        for {set xi 0} { $xi < $dim } {incr xi} {
             lset rvec $j [expr [lindex $rvec $j] + [lindex $vec $xi] * [lindex $hvec $j $xi]]
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

#############################################################
#
# Prepare Connection to VMD
#
#############################################################

proc prepare_vmd_connection { {filename "vmd"} {wait "0"} {start "1" } } {
    writepsf "$filename.psf"
    writepdb "$filename.pdb"
    for {set port 10000} { $port < 65000 } { incr port } {
        catch {imd connect $port} res
        if {$res == ""} break
    }
    set HOSTNAME [exec hostname]
    set vmdout_file [open "vmd_start.script" "w"]
    puts $vmdout_file "mol load psf $filename.psf pdb $filename.pdb"
    puts $vmdout_file "rotate stop"
    puts $vmdout_file "mol modstyle 0 0 CPK 1.000000 0.300000 8.000000 6.000000"
    puts $vmdout_file "mol modcolor 0 0 SegName"
    puts $vmdout_file "imd connect $HOSTNAME $port"
     close $vmdout_file
    if { $start == 0 } {
        puts "Start VMD in the same directory on the machine you with :"
        puts "vmd -e vmd_start.script &"
        imd listen $wait
    } else {
        exec vmd -e vmd_start.script &
    }
}

# System identification: 
set name  "shaved-PPV"
set ident "temp"
set mypi          3.141592653589793


#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#############################################################
#                                                           #
#  Template for Shaved-Bundle Simulations                   #
#                                                           #
#                                                           #
#  Created:       23.01.2003 by Sayar                       #
#  Last modified: 23.01.2003 by Sayar                       #
#                                                           #
#############################################################

puts " "
puts "======================================================="
puts "=                shaved-bundle_temp.tcl               ="
puts "= This script tests the following features in Espresso="
puts "= -- lj-cos, lj, coulomb (unscreened, non-periodic),  ="
puts "=    bond and angle potentials                        ="
puts "= -- constraints (wall, cylinder,sphere)              ="
puts "======================================================="
puts " "

# System parameters
#############################################################

set bundle_length 1
set n_poly 7
set real_p_length 5

set p_length [expr $real_p_length * $bundle_length ]
set b_length 1.00
set d_max 1.12
#set d_space {[expr $d_max] [expr $d_max] 1}
set d_space { 0 0 0}
lset d_space 0 [expr $d_max]
lset d_space 1 [expr $d_max]
lset d_space 2 1

set p_val 1
set c_ion_val 1.0
set density 0.0001 


# Interaction parameters
#############################################################

set ljr_cut       1.12246204831
set ljr_eps       1.0


set ljcosa_cut       1.5

set ljcosa_eps       4.0
set ljcosa_eps_end   4.0

set ljcosa_sig       1.0
set ljcosa_rmin 1.12246204831


set fene_r        2.0
set fene_k        7.0

set bend_k        10.0

set bjerrum       0.0
set cion_d      1.0

set accuracy      1.0e-4

# Integration parameters
#############################################################

set time_step    0.005
set skin         0.5

set warm_steps   1000
set warm_n_times 1

set min_dist     1.1

set ci_steps     1000
set ci_n_times   1

set poly_steps     1000
set poly_n_times   1

set int_steps    1000
set int_n_times  1

# Other parameters
#############################################################
set tcl_precision 6
set mypi          3.141592653589793
set vmd_output    "no"
set vmd_wait      0

#############################################################
#  Setup System                                             #
#############################################################

#set density       0.0001
set n_part [expr $n_poly * $p_length * (1. + 1./$c_ion_val) ]
set volume [expr $n_part/$density]
set sphere_rad [expr pow((3.0*$volume)/(4.0*$mypi),1.0/3.0)]
set  box_l       [expr 4.0*$sphere_rad + 6.0*$skin]
set center      [expr $box_l/2.0]
set shift_vec {0 0 0}
lset shift_vec 0 [expr $box_l/2.0] 
lset shift_vec 1 [expr $box_l/2.0] 
lset shift_vec 2 [expr $box_l/2.0- ($p_length-1.)/2.*$b_length]

setmd box_l     $box_l $box_l $box_l
setmd periodic  0 0 0
setmd time_step $time_step
setmd skin      $skin
thermostat langevin 1.0 1.0
cellsystem NSQUARE

# Interaction setup
#############################################################

# type 0 charged backbone, type 1 backbone chain_end, type 2 counterion
set max_inter_type 3

# repulsive LJ for all
set ljr_shift  0.25
inter 0 0 lennard-jones $ljr_eps 1.0 $ljr_cut $ljr_shift 0
inter 0 1 lennard-jones $ljr_eps 1.0 $ljr_cut $ljr_shift 0
inter 1 1 lennard-jones $ljr_eps 1.0 $ljr_cut $ljr_shift 0

inter 0 2 lennard-jones $ljr_eps $cion_d $ljr_cut $ljr_shift 0
inter 1 2 lennard-jones $ljr_eps $cion_d $ljr_cut $ljr_shift 0

inter 2 2 lennard-jones $ljr_eps $cion_d $ljr_cut $ljr_shift 0


# FENE 
inter 0 fene $fene_k $fene_r

# Stiffness
inter 1 angle $bend_k

# Constraint Sphere
constraint sphere center $center $center $center radius $sphere_rad type 3 
inter 0 3 lennard-jones 1.0 1.0 1.12246 0.25 0 0
inter 1 3 lennard-jones 1.0 1.0 1.12246 0.25 0 0
inter 2 3 lennard-jones 1.0 1.0 1.12246 0.25 0 0


#############################################################
#  Create a bundle of n_poly chains                         #
#############################################################

set part_id 0
set input_file [open "pe_micelle.data" "r"]
for {set xi 0} { $xi < [expr $n_poly] } {incr xi} {
    gets $input_file vec
    set rvec [hexconvert $vec $shift_vec $d_space]
    set part_id [create_chain $part_id $rvec $p_length $b_length]
}
set d_space {2. 2. 1}
for {set xi 0} { $xi < [expr $n_poly] } {incr xi} {
    gets $input_file vec
    set rvec [hexconvert $vec $shift_vec $d_space]
    set part_id [create_chain $part_id $rvec [expr $p_length/$c_ion_val] $b_length]
}
close $input_file


#fix the particles on the chains and set up lj types and add stiffness and put the bonds
for {set xi 0} { $xi < [expr $n_poly * $p_length ]} {incr xi} {
    if { ([expr $xi % $real_p_length ]  == 0) || ([expr $xi % $real_p_length ]  == [expr $real_p_length -1 ])} {
        part $xi type 1 q $p_val fix
    } else {
        part $xi type 0 q $p_val fix
    }
	if { [expr $xi % $real_p_length ]  != 0 } {
        part [expr $xi-1] bond 0 $xi
        if { [expr $xi % $real_p_length]  != 1 } {
            part [expr $xi-1] bond 1 [expr $xi-2] [expr $xi]
        }
	}
}

# particle numbers
set n_mono [expr $n_poly * $p_length]
set n_ci   [expr $n_mono/$c_ion_val]
set n_part [expr $n_mono + $n_ci ]

# Counterions 
for {set xi 0} { $xi < $n_ci } {incr xi} {
    part [expr $n_mono + $xi]  type 2 q -$c_ion_val
}

#now add a cylindrical constraint at the center of mass of the particles

# calculate center of mass
set temp0 0.
set temp1 0.
set temp2 0.
set bundle_cm { 0. 0. 0. }

for {set xi 0} { $xi < [expr $n_poly * $p_length ] } {incr xi} {
    for {set j 0} { $j < 3 } {incr j} {
        lset bundle_cm $j [expr [lindex $bundle_cm $j]  + [lindex [part $xi] [expr ($j + 2)] ]]
    }
}

for {set j 0} { $j < 3 } {incr j} {
    lset bundle_cm $j [expr [lindex $bundle_cm $j] / ($n_poly * $p_length )]
}


set cylinder_skin 1.20
set cylinder_cutoff 1.12246204831
set cylinder_rad [expr ((sqrt ($n_poly )) * $ljcosa_rmin * .5 * $d_max)  ]
set cylinder_half_length [expr  ($p_length - 1.) / 2. * $b_length ]

constraint cylinder center [lindex $bundle_cm 0] [lindex $bundle_cm 1] [lindex $bundle_cm 2] axis 0 0 1 radius [expr $cylinder_rad+ $cylinder_cutoff * $cylinder_skin] length [expr $cylinder_half_length+ $cylinder_cutoff * $cylinder_skin] direction -1 type 4; inter 0 4 lennard-jones 1.0 1.0 1.12246204831 0.25 0 0; inter 1 4 lennard-jones 1.0 1.0 1.12246204831 0.25 0 0

#add a cylinder for the counter ion equilibreation. 
#puts "Add a cylinder for the counter ion equilibreation."
set skin_ci_fac 5.0

#constraint cylinder center [lindex $bundle_cm 0] [lindex $bundle_cm 1] [lindex $bundle_cm 2] axis 0 0 1 radius [expr $cylinder_rad +$cylinder_cutoff * $skin_ci_fac] length [expr $cylinder_half_length +$cylinder_cutoff * $skin_ci_fac] direction -1 type 5; inter 2 5 lennard-jones 1.0 1.0 1.12246204831 0.25 0 0

constraint cylinder center [lindex $bundle_cm 0] [lindex $bundle_cm 1] [lindex $bundle_cm 2] axis 0 0 1 radius [expr $cylinder_rad +$cylinder_cutoff * $skin_ci_fac] length 500 direction -1 type 5; inter 2 5 lennard-jones 1.0 1.0 1.12246204831 0.25 0 0

constraint wall normal 0 0 -1 dist [expr -[lindex $bundle_cm 2]-[expr $cylinder_half_length +$cylinder_cutoff * $skin_ci_fac]] type 8; inter 2 8 lennard-jones 1.0 1.0 1.12246204831 0.25 0 0
constraint wall normal 0 0 1 dist [expr [lindex $bundle_cm 2]-[expr $cylinder_half_length +$cylinder_cutoff * $skin_ci_fac]] type 9; inter 2 9 lennard-jones 1.0 1.0 1.12246204831 0.25 0 0

#add a cylinder to stop c_ions penetrating into the bundle. 
#puts "Add a cylinder to stop c_ions penetrating into the bundle."
set skin_ci_fac -0.5

constraint cylinder center [lindex $bundle_cm 0] [lindex $bundle_cm 1] [lindex $bundle_cm 2] axis 0 0 1 radius [expr $cylinder_rad +$cylinder_cutoff * $skin_ci_fac] length [expr $cylinder_half_length +$cylinder_cutoff * $skin_ci_fac] direction 1 type 6; inter 2 6 lennard-jones 1.0 1.0 1.12246204831 0.25 0 0

#puts "Constraints:\n[constraint]"


#for {set i 0} { $i < $part_id } {incr i} {puts "[part $i]"}

# Status report
#puts "Simulate $n_part particles in a spherical cell with radius $sphere_rad"
#puts "$n_poly PE backbones of length $p_length and charge distance $b_length"
#puts "Each monomer has $p_val charge"
#puts "neutralized by $n_ci counterions with valency $c_ion_val."
#puts "Constraints:\n[constraint]"


#############################################################
#  VMD connection                                           #
#############################################################
if { $vmd_output=="yes" } {
    prepare_vmd_connection "$name$ident" $vmd_wait 1
}




#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#############################################################
#  Warmup Integration                                       #
#############################################################

#puts "\nStart warmup integration: max $warm_n_times times $warm_steps steps\n"
set act_min_dist [analyze mindist]
set cap 20  
inter forcecap $cap

for {set xi 0} { $xi < $warm_n_times  && $act_min_dist < $min_dist } { incr xi} {
    puts -nonewline "Run $xi at time=[setmd time] min_dist=$act_min_dist\r"
    flush stdout
    integrate $warm_steps
    if { $vmd_output=="yes" } {catch { imd positions }}
    set act_min_dist [analyze mindist]
    set cap [expr $cap+10]
    inter forcecap $cap  
}

#puts "Warmup done: Minimal distance [analyze mindist]\n"

#############################################################
#  Counterion Equilibration                                #
#############################################################

#puts "Start Counterion Equilibration: $ci_n_times times $ci_steps steps\n"
# remove LJ cap
inter forcecap 0
# Start Coulomb interaction
inter coulomb $bjerrum dh 0.0 [expr 2.0*$sphere_rad]

#puts "Interactions:\n [inter]"

for {set xi 0} { $xi < $ci_n_times  } { incr xi} {
    puts -nonewline "Run $xi at time=[setmd time]\r"
    flush stdout
    integrate $ci_steps
    if { $vmd_output=="yes" } { imd positions }
}

#puts "Counterion Equilibration done."

#release the counter ions out of the cylinder 
inter 2 5 lennard-jones 0.0 1.0 0.0 0.25 0 0

#############################################################
#      Polymer equilibration within the cylinder constraint #
#############################################################
# free backbone monomers
for {set xi 0} { $xi < $n_mono } {incr xi} {
    part $xi v 0 0 0 f 0 0 0 unfix 
}

#change the interaction type for backbone to lj+cos
#first set the repulsive lj cutoff to zero
set ljr_eps 0
set ljr_cut 0
inter 0 0 lennard-jones $ljr_eps 1.0 $ljr_cut 0.25 0
inter 0 1 lennard-jones $ljr_eps 1.0 $ljr_cut 0.25 0
inter 1 1 lennard-jones $ljr_eps 1.0 $ljr_cut 0.25 0
#next put the attractive lj+cos
inter 0 0 lj-cos $ljcosa_eps $ljcosa_sig $ljcosa_cut 0.0
inter 0 1 lj-cos $ljcosa_eps $ljcosa_sig $ljcosa_cut 0.0
inter 1 1 lj-cos $ljcosa_eps_end $ljcosa_sig $ljcosa_cut 0.0

for {set xi 0} { $xi < $poly_n_times } { incr xi} {
    puts -nonewline "Run $xi at time=[setmd time]\r"
    flush stdout
    integrate $poly_steps
    if { $vmd_output=="yes" } { imd positions }
}

#puts "Polymer Equilibration done."

#Relax the wall constraints on the polymers 
inter 0 4 lennard-jones 0.0 1.0 0.0 0.25 0 0
inter 1 4 lennard-jones 0.0 1.0 0.0 0.25 0 0

#puts "Constraints= [constraint]"

#puts "Interactions= [inter]"

#############################################################
#      Integration                                          #
#############################################################



#reset the forces and velovities for all the particles.
kill_particle_motion

analyze set chains 0 $n_poly $p_length
#puts "Interactions:\n [inter]"
#puts "Constraints= [constraint]\n\n"

#puts "\nStart integration: $int_n_times times $int_steps steps"
#set energyfile [open "energy.txt" "w"]

for {set xi 0} { $xi < $int_n_times } { incr xi} {
    #puts  "Run $xi at time=[setmd time], re = [analyze re], Energy = [analyze energy]"
    integrate $int_steps
    if { $vmd_output=="yes" } { imd positions }

    #set out [open "|gzip -c - > $name.[format %05d $xi].gz" "w"]
    #blockfile $out write particles "id pos type q v f" all
    #blockfile $out write bonds all
    #close $out
}

#puts "\nIntegration done."

puts  "Energy = [analyze energy]"

puts "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nCOMPARE ENERGY VALUES ABOVE WITH THE ONES BELOW\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n Energy = { energy -216.909 } { kinetic 92.9468 } { 0 FENE 141.601 } { 1 angle 6.56779 } { 0  0 nonbonded -248.561 } { 0  1 nonbonded -139.074 } { 0  2 nonbonded 0.0 } { 0  3 nonbonded 0.0 } { 1  1 nonbonded -70.3895 } { 1  2 nonbonded 0.0 } { 1  3 nonbonded 0.0 } { 2  2 nonbonded 0.0 } { 2  3 nonbonded 0.0 } { 2  6 nonbonded 0.0 } { 2  8 nonbonded 0.0 } { 2  9 nonbonded 0.0 }\n"

puts "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\nIF THE VALUES ARE CORRECT, IT IS YOUR LUCKY DAY !!!\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n"


#analyze append; 
#checkpoint_set "checkpoint.gz" 

#writepsf "output.psf"
#writepdb "output.pdb"
#############################################################

puts "\nFinished"
exit 0
