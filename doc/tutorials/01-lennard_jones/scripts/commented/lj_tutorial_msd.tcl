# Copyright (C) 2010,2011,2012 The ESPResSo project
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

#############################################################
#                                                           #
#  Calculate the velocity autocorrelation function          # 
#  using the Espresso routines                              #
#                                                           #
#############################################################

source lj_functions.tcl
cellsystem domain_decomposition -no_verlet_list

set n_part 108

set lj1_eps     1.0
set lj1_sig     1.0
set lj1_cut     2.5
set lj1_shift   [expr -(pow(1.0/$lj1_cut,12)-pow(1.0/$lj1_cut,6))]
set lj1_offset  0.0

thermostat off

set eq_tstep 0.0001
set tstep 0.001
set skin 0.1
setmd skin $skin
set target_temperature 0.728

set warm_steps   100
set warm_n_times 2000
set min_dist     0.87

set sampling_interval 10
set equilibration_interval 1000
set sampling_iterations 10000
set equilibration_iterations 200

set tcl_precision 8
t_random seed 12345

set density 0.8442
set box_length [expr pow($n_part/$density,1.0/3.0)+2*$skin]
setmd box $box_length $box_length $box_length

for {set i 0} {$i < $n_part} {incr i} {
   set pos_x [expr [t_random]*$box_length]
   set pos_y [expr [t_random]*$box_length]
   set pos_z [expr [t_random]*$box_length]
   part $i pos $pos_x $pos_y $pos_z q 0.0 type 0
}

inter 0 0 lennard-jones $lj1_eps $lj1_sig $lj1_cut $lj1_shift $lj1_offset

set act_min_dist [analyze mindist]ity
set cap 1.0
inter forcecap $cap

setmd time_step $tstep
set i 0
while { $i < $warm_n_times && $act_min_dist < $min_dist } {
    integrate $warm_steps

    set act_min_dist [analyze mindist]
    set cap [expr $cap+1.0]
    inter forcecap $cap
    incr i
}

inter forcecap 0
setmd time_step $eq_tstep

for { set i 0 } { $i < $equilibration_iterations } { incr i } {
   integrate $equilibration_interval
   rescale_velocities $target_temperature $n_part
}

setmd time_step $tstep
set tmax 25.0
set delt [expr $tstep*$sampling_interval]
set pos_id [observable new particle_positions type 0]
set msd [correlation new obs $pos_id corr_operation square_distance_componentwise dt $delt tau_max $tmax]

for {set i 0} { $i < $sampling_iterations } { incr i} {
     integrate $sampling_interval
     correlation $msd update
}

correlation $msd finalize
correlation $msd write_to_file "data/msd_internal.dat"

exit
