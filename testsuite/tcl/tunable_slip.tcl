# Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#   
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
#                                                           #
# Basic tests of the tunable-slip boundary interactions     #
#                                                           #
# 1) check constraints                                      #
# 2) check implementation with DPD                          #
#                                                           #
# Reference:                                                #
# J.Smiatek, M.P. Allen, F. Schmid:                         #
# "Tunable-slip boundaries for coarse-grained simulations   #
# of fluid flow", Europ. Phys. J. E 26, 115 (2008)          #
#                                                           #
#############################################################

source "tests_common.tcl"

require_feature "TUNABLE_SLIP"
require_feature "DPD"
require_feature "CONSTRAINTS"
require_feature "EXTERNAL_FORCES"
require_feature "LENNARD_JONES"

puts "----------------------------------------"
puts "- Testcase tunable_slip.tcl running on [format %02d [setmd n_nodes]] nodes  -"
puts "----------------------------------------"

# System parameters 
set box_l                10.0
set box_x                $box_l
set box_y                $box_l
set box_z                $box_l

# skin depth, not used here
setmd skin               0.0

# box length
setmd box_l              $box_l $box_l $box_l

# volume of system
set volume               [expr $box_x*$box_y*$box_z]

# periodic boundary conditions
setmd periodic           1 1 1

################################################
#         Solvent parameters                   #
################################################

# type number for solvent molecules
set solvent_id           0

# density for solvent molecules
set density_s            3.0

# number of solvent particles
set n_solvent            [expr int (ceil($density_s*$volume))]

# total number of particles
set n_part               $n_solvent

################################################
#         Interaction parameters               #
################################################

# Lennard-Jones ################################

# epsilom
set lj_eps               1.0

# sigma
set lj_sig               1.0

# cutoff
set lj_cut               1.12246

# shift
set lj_shift             0.0

# offset
set lj_off               0.0

################################################
#            Thermostat                        #
################################################

# temperature
set temp                 1.0

#Tunable Slip Boundaries
set gamma_L              1.0

# cutoff has to be larger than LJ-interaction range
set r_cut_L              2.0

# Extrem accelaration of the flow
set vx                   10.0

set vy                   0.0

set vz                   0.0

# DPD thermostat ###############################

set dpd_temperature            1.0

# DPD friction coefficient
set dpd_gamma                  5.0

# DPD-r_cutoff radius
set dpd_r_cut                  1.0

# length of time step
setmd time_step                 0.01

set timestep                    0.01

# Length of integration scheme
set int_steps                   100

# Length of integration run
set int_loops                   3

# precision of data output
set tcl_precision     6

# initialize random number generator
set t_random_seed     [expr int(rand()*99999)^[pid]]

# external force
set f_x    1.0
set f_y    0.0
set f_z    0.0

cellsystem domain_decomposition -no_verlet_list

# Constraints
set wall_left_id 1
set wall_right_id 2

constraint wall normal 0 0  1 dist   0 type $wall_left_id
constraint wall normal 0 0 -1 dist -10 type $wall_right_id

# DPD-Thermostat
thermostat dpd $dpd_temperature $dpd_gamma $dpd_r_cut

# Solvent particles
for {set i 0} { $i < $n_solvent} { incr i } {
# Placing particles inside the constraints
    set posx [expr 1.0+8.0*[t_random]]
    set posy [expr 1.0+8.0*[t_random]]
    set posz [expr 1.0+8.0*[t_random]]

    set vx   [expr 0.1*[t_random]]
    set vy   [expr 0.1*[t_random]]
    set vz   [expr 0.1*[t_random]]
    
    part $i pos $posx $posy $posz type $solvent_id v $vx $vy $vz ext_force $f_x $f_y $f_z
}
galilei_transform

# Interactions 
inter $wall_right_id $solvent_id lennard-jones $lj_eps $lj_sig $lj_cut $lj_shift $lj_off 
inter $wall_left_id $solvent_id lennard-jones $lj_eps $lj_sig $lj_cut $lj_shift $lj_off 
inter $wall_left_id $solvent_id tunable_slip $temp $gamma_L $r_cut_L $timestep $vx $vy $vz
inter $wall_right_id $solvent_id tunable_slip $temp $gamma_L $r_cut_L $timestep $vx $vy $vz

############################################
#      Procedures                          #
############################################

proc measure_kinetic_energy {} {
    global n_solvent
    global energy
    global n_solvent 
    global energy
    global E_ref

    set energy [analyze energy kinetic]
    set E_ref [expr $n_solvent*1.5]
    if {$energy <= $E_ref} {
	error "Tunable-slip layer does not work ..."
    }
}

################### Integration #############################

puts "cells = [setmd cell_grid]"
puts "max_range = [setmd max_range]"
puts "n_particles = [setmd n_part]"

if { [catch {

    for {set step 0} {$step < $int_loops} {incr step} {
	puts "step $step"
	integrate $int_steps
    }

    measure_kinetic_energy

} res ] } {
    error_exit $res
}

exit 0
