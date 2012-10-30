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
#  Sample System  Lennard Jones Liquid                      #
#                                                           #
#    LJ System Studied in the Book                          # 
#    Understanding Molecular Simulations                    # 
#    Frankel & Smith , 2nd Edition                          # 
#    Chapter 4, Case Study 4                                # 
#                                                           # 
#############################################################

# Source (call) functions to be used
source lj_functions.tcl

# Shows the incompiled features of espresso
puts " "
puts "========================================================"
puts "=                lj_liquid_tutorial.tcl                ="
puts "========================================================"
puts " "
puts "Espresso Code Base : \n[code_info]\n"
puts "========================================================"
puts " "

# Use a cell list, but not a Verlet list
cellsystem domain_decomposition -no_verlet_list

#############################################################
#                                                           #
#  Parameters                                               #
#                                                           #
#############################################################

# System parameters
#############################################################

# 108 particles
set n_part 108

# Interaction parameters
#############################################################

# The energy, particle 'size', cut-off radius, potential 
# shift, and potential offset. The shift is such that the 
# potential is zero at the cut-off radius
set lj1_eps     1.0
set lj1_sig     1.0
set lj1_cut     2.5
set lj1_shift   [expr -(pow(1.0/$lj1_cut,12)-pow(1.0/$lj1_cut,6))]
set lj1_offset  0.0

# Integration parameters
#############################################################

# Simulation in NVE Ensemble
thermostat off

# Define the time step during equilibration, the regular time 
# step, the size of the cell-list skin + put this value in
# ESPResSo, and the target temperature
set eq_tstep 0.0001
set tstep 0.001
set skin 0.1
setmd skin $skin
set target_temperature 0.728

# Define the number of warmup steps and times. Require that
# the particles have at least the distance min_dist after
set warm_steps   100
set warm_n_times 2000
set min_dist     0.87

# Production/Equilibration values for the integration length
set sampling_interval 10
set equilibration_interval 1000
set sampling_iterations 10000
set equilibration_iterations 200

# Other parameters
#############################################################

# precision for the puts output
set tcl_precision 8

# setting a seed for the random number generator
t_random seed 12345

#############################################################
#                                                           #
#  Setup System                                             #
#                                                           #
#############################################################

# Box setup
#############################################################

set density 0.8442

# determine the box sizes from the density and set
set box_length [expr pow($n_part/$density,1.0/3.0)+2*$skin]
puts "density = $density, box_length = $box_length"
setmd box $box_length $box_length $box_length

# Particle setup
#############################################################

# set the particle idenifier, position, charge and type 
for {set i 0} {$i < $n_part} {incr i} {
  set pos_x [expr [t_random]*$box_length]
  set pos_y [expr [t_random]*$box_length]
  set pos_z [expr [t_random]*$box_length]
  part $i pos $pos_x $pos_y $pos_z q 0.0 type 0
}

# Interaction setup
#############################################################

inter 0 0 lennard-jones $lj1_eps $lj1_sig $lj1_cut $lj1_shift $lj1_offset

set act_min_dist [analyze mindist]
puts "Start with minimal distance $act_min_dist"

#############################################################
#                                                           #
#  Warmup Integration                                       #
#                                                           #
#############################################################

puts "\nStart warmup integration:"
puts "At maximum $warm_n_times times $warm_steps steps"
puts "Stop if minimal distance is larger than $min_dist"

# set LJ cap, i.e., the forcecap which prevents (partially)
# overlapping particles from being propelled at great velocity
set cap 1.0
inter forcecap $cap

# Warmup integration loop, in which the capping is slowly increased
setmd time_step $tstep
set i 0
while { $i < $warm_n_times && $act_min_dist < $min_dist } {
  integrate $warm_steps

  # Warmup criterion
  set act_min_dist [analyze mindist]
  puts -nonewline "run $i at time = [setmd time] (LJ cap = $cap) min dist = $act_min_dist\r"
  flush stdout

  # Increase LJ cap
  set cap [expr $cap+1.0]
  inter forcecap $cap
  incr i
}
# It is assumed that the number of $warm_n_times*$warm_steps of 
# integration is sufficient to reach the desired minimum distance.
# This is true for the parameters we chose here

puts "\nWarmup run finished"

# disable force capping when done
inter forcecap 0

#############################################################
#                                                           #
#  Equilibration and Production                             #
#                                                           #
#############################################################

# set the time step for the equilibration part
setmd time_step $eq_tstep

puts "\nStart equilibration integration:"
for { set i 0 } { $i < $equilibration_iterations } { incr i } {
  integrate $equilibration_interval
  rescale_velocities $target_temperature $n_part
  puts -nonewline "equilibration run $i at time = [setmd time]\r"
  flush stdout
}
puts "\nEquilibration run finished"

# set output file for information on the simulation
set blockfile [open "data/sim_info.dat" "w"]

# set output file for the energy and make a header for it
set en [open "data/energy.dat" "w"]
puts $en "# Iteration Pressure Kinetic Potential Temperature Total"

# set the time step for the production part
setmd time_step $tstep

puts "\nStart production integration:"
# start sampling loop
set nrep 0
for {set i 0} { $i < $sampling_iterations } { incr i } {
  # counter to keep track of the total number of cycles
  incr nrep 

  # perform the MD integration
  integrate $sampling_interval

  # take all the particle coordinates, velocities etc,
  # and put them in the blockfile format and use the
  # save_sim routine from the lf_functions.tcl file 
  # to save the data to "data/sim_info.dat"
  save_sim $blockfile "id pos v f q type" "all"

  # obtain all relevant physical quantities and write to
  # the "data/energy.dat" file. The energies are per particle

  set energies [analyze energy]
  set pressure [analyze pressure total]
  set total [expr [lindex $energies 0 1]/$n_part]
  set kinetic [expr [lindex $energies 1 1]/$n_part]
  set potential [expr [lindex $energies 2 3]/$n_part]
  set temperature [expr $kinetic/(1.5)]
  puts $en "$i $pressure $kinetic $potential $temperature $total"

  # maintain lists for the physical quantities during the 
  # production runs
  lappend apressure $pressure
  lappend akinetic $kinetic
  lappend apotential $potential
  lappend atemperature $temperature
  lappend atotal $total

  # show progress bar
  puts -nonewline "integration step $i / $sampling_iterations\r"
  flush stdout
}
puts "\nProduction run finished \n"

# writing to file is done
close $en
close $blockfile

#############################################################
#                                                           #
#  Statistics and Data                                      #
#                                                           #
#############################################################

# now calculate mean, the statistical error, and the error in
# the error using the 'uwerr' command, which employs the Wolff
# method, for all relevant physical quantities

puts "Reporting Energies and Temperature"

set error [uwerr $atotal $nrep 1 ]
set value [lindex $error 0]
set verror [lindex $error 1]
puts "  Total Energy: $value  $verror"

set error [uwerr $akinetic $nrep 1 ]
set value [lindex $error 0]
set verror [lindex $error 1]
puts "  Kinetic Energy: $value  $verror"

set error [uwerr $apotential $nrep 1 ]
set value [lindex $error 0]
set verror [lindex $error 1]
puts "  Potential Energy: $value  $verror"

set error [uwerr $atemperature $nrep 1 ]
set value [lindex $error 0]
set verror [lindex $error 1]
puts "  Temperature : $value  $verror"

set error [uwerr $apressure $nrep 1 ]
set value [lindex $error 0]
set verror [lindex $error 1]
puts "  Pressure : $value  $verror"

# Congratulations, you have now successfully executed
# your first Espresso simulation
puts "========================================================"
exit
