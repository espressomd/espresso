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

source lj_functions.tcl

puts " "
puts "========================================================"
puts "=                lj_liquid_tutorial.tcl                ="
puts "========================================================"
puts " "
puts "Espresso Code Base : \n[code_info]\n"
puts "========================================================"
puts " "

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
puts "density = $density, box_length = $box_length"
setmd box $box_length $box_length $box_length

for {set i 0} {$i < $n_part} {incr i} {
  set pos_x [expr [t_random]*$box_length]
  set pos_y [expr [t_random]*$box_length]
  set pos_z [expr [t_random]*$box_length]
  part $i pos $pos_x $pos_y $pos_z q 0.0 type 0
}

inter 0 0 lennard-jones $lj1_eps $lj1_sig $lj1_cut $lj1_shift $lj1_offset

set act_min_dist [analyze mindist]
puts "Start with minimal distance $act_min_dist"

puts "\nStart warmup integration:"
puts "At maximum $warm_n_times times $warm_steps steps"
puts "Stop if minimal distance is larger than $min_dist"

set cap 1.0
inter forcecap $cap

setmd time_step $tstep
set i 0
while { $i < $warm_n_times && $act_min_dist < $min_dist } {
  integrate $warm_steps

  set act_min_dist [analyze mindist]
  puts -nonewline "run $i at time = [setmd time] (LJ cap = $cap) min dist = $act_min_dist\r"
  flush stdout

  set cap [expr $cap+1.0]
  inter forcecap $cap
  incr i
}

puts "\nWarmup run finished"

inter forcecap 0

setmd time_step $eq_tstep

puts "\nStart equilibration integration:"
for { set i 0 } { $i < $equilibration_iterations } { incr i } {
  integrate $equilibration_interval
  rescale_velocities $target_temperature $n_part
  puts -nonewline "equilibration run $i at time = [setmd time]\r"
  flush stdout
}
puts "\nEquilibration run finished"

set blockfile [open "data/sim_info.dat" "w"]

set en [open "data/energy.dat" "w"]
puts $en "# Iteration Pressure Kinetic Potential Temperature Total"

setmd time_step $tstep

puts "\nStart production integration:"

set nrep 0
for {set i 0} { $i < $sampling_iterations } { incr i } {

  incr nrep 

  integrate $sampling_interval

  save_sim $blockfile "id pos v f q type" "all"

  set energies [analyze energy]
  set pressure [analyze pressure total]
  set total [expr [lindex $energies 0 1]/$n_part]
  set kinetic [expr [lindex $energies 1 1]/$n_part]
  set potential [expr [lindex $energies 2 3]/$n_part]
  set temperature [expr $kinetic/(1.5)]
  puts $en "$i $pressure $kinetic $potential $temperature $total"

  lappend apressure $pressure
  lappend akinetic $kinetic
  lappend apotential $potential
  lappend atemperature $temperature
  lappend atotal $total

  puts -nonewline "integration step $i / $sampling_iterations\r"
  flush stdout
}
puts "\nProduction run finished \n"

close $en
close $blockfile

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

puts "========================================================"
exit
