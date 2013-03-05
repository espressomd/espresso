# Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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


puts " "
puts "======================================================="
puts "=                        lj_liquid_tutorial.tcl       ="
puts "======================================================="
puts " "
puts "Espresso Code Base : \n[code_info]\n"
puts " "
puts " "



cellsystem domain_decomposition  -no_verlet_list
#############################################################
#  Parameters                                               #
#############################################################
# System identification:
set name  "lj_liquid"
set ident "_s1"
# System parameters
#############################################################

# 108 particles
set n_part  108

# Interaction parameters
#############################################################

set lj1_eps     1.0
set lj1_sig     1.0
set lj1_cut     2.5
# set lj1_shift   [expr -4.0*$lj1_eps*(pow(1.0/$lj1_cut,12)-pow(1.0/$lj1_cut,6))]
set lj1_shift   [expr -(pow(1.0/$lj1_cut,12)-pow(1.0/$lj1_cut,6))]
set lj1_offset  0.0

# Integration parameters
#############################################################
thermostat off ;# Simulation in NVE Ensemble

setmd time_step 0.001
set eq_tstep 0.0001
set tstep   0.001
set skin 0.1
setmd skin  $skin
set target_temperature 0.728

set warm_steps   100
set warm_n_times 2000
# do the warmup until the particles have at least the distance min__dist
set min_dist     0.87

# integration
set sampling_interval    10
set equilibration_interval 1000

set sampling_iterations 10000
set equilibration_iterations 200

# Other parameters
#############################################################
set tcl_precision 8

#setting a seed for the random number generator
expr srand([pid])

#############################################################
#  Setup System                                             #
#############################################################


# Particle setup
#############################################################

set density 0.8442
set box_length [expr pow($n_part/$density,1.0/3.0)+2*$skin]
puts " density = $density box_length = $box_length"
setmd box $box_length $box_length $box_length
for {set i 0} {$i < $n_part} {incr i} {
   set pos_x [expr rand()*$box_length]
   set pos_y [expr rand()*$box_length]
   set pos_z [expr rand()*$box_length]
   part $i pos $pos_x $pos_y $pos_z q 0.0 type 0
}

writepdb data/config.pdb

# Interaction setup
#############################################################
inter 0 0 lennard-jones $lj1_eps $lj1_sig $lj1_cut $lj1_shift $lj1_offset

set act_min_dist [analyze mindist]
puts "Start with minimal distance $act_min_dist"

#############################################################
#  Warmup Integration                                       #
#############################################################


# open Observable file

puts "\nStart warmup integration:"
puts "At maximum $warm_n_times times $warm_steps steps"
puts "Stop if minimal distance is larger than $min_dist"

# set LJ cap
set cap 1.0
inter forcecap $cap

# Warmup Integration Loop
set i 0
while { $i < $warm_n_times && $act_min_dist < $min_dist } {

    integrate $warm_steps

    # Warmup criterion
    set act_min_dist [analyze mindist]
    puts -nonewline "run $i at time=[setmd time] (LJ cap=$cap) min dist = $act_min_dist\r"
    flush stdout

#   write observables

#   Increase LJ cap
    set cap [expr $cap+1.0]
    inter forcecap $cap
    incr i
}

inter forcecap 0

# 
puts "\n Warm up finished \n"
setmd time_step $eq_tstep
for { set i 0 } { $i < $equilibration_iterations } { incr i } {
   integrate $equilibration_interval
   set energies [analyze energy]
   rescale_velocities  $target_temperature [setmd n_part]
   puts -nonewline "eq run $i at time=[setmd time] \r"
   # puts "temp = [expr [lindex $energies 1 1]/(1.5*[setmd n_part])]"
}
  puts "sampling "
setmd time_step $tstep
set nrep 0

set en [open "data/energy.dat" "w"]
set blockfile [open "data/sim_info.dat" "w"]

puts $en "#"
puts $en "#"
puts $en "#  Pressure    Kinetic    Potential  Temperature "
puts $en "# "
for {set i 0} { $i < $sampling_iterations } { incr i} {
      incr nrep 
     integrate $sampling_interval
     save_sim $blockfile "id pos v f q type" "all"
     set energies [analyze energy]
     set pressure [analyze pressure total]
   # puts "pressure=$pressure "
   # puts "energies: $energies"
     set total [expr [lindex $energies 0 1]/$n_part]
     set kinetic [expr [lindex $energies 1 1]/$n_part]
     set potential [expr [lindex $energies 2 3]/$n_part]
     set kinetic_temperature [expr [lindex $energies 1 1]/(1.5*[setmd n_part])]
     set temperature $kinetic_temperature
     puts $en "$i   $pressure  $kinetic  $potential  $temperature $otal"
     lappend apressure $pressure
     lappend akinetic $kinetic
     lappend apotential $potential
     lappend atemperature $temperature
     lappend atotal $total
     puts -nonewline "integration step   $i / $sampling_iterations \r"
     # puts -nonewline "total = [expr $total/$n_part] kinetic= [expr $kinetic/$n_part] potential= [expr $potential/$n_part]  temperature = $temperature\n"
}
close $en
  puts "--Reporting Energies and Temperature"
  set error [uwerr $atotal $nrep 1 ]
  set value [lindex $error 0]
  set verror [lindex $error 1]
  puts "     Total Energy: $value  $verror"
  set error [uwerr $akinetic $nrep 1 ]
  set value [lindex $error 0]
  set verror [lindex $error 1]
  puts "     Kinetic Energy: $value  $verror"
  set error [uwerr $apotential $nrep 1 ]
  set value [lindex $error 0]
  set verror [lindex $error 1]
  puts "     Potential Energy: $value  $verror"
  set error [uwerr $atemperature $nrep 1 ]
  set value [lindex $error 0]
  set verror [lindex $error 1]
  puts "     Temperature : $value  $verror"
  set error [uwerr $apressure $nrep 1 ]
  set value [lindex $error 0]
  set verror [lindex $error 1]
  puts "     Pressure : $value  $verror"
exit
