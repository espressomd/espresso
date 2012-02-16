#
# Copyright (C) 2010,2012 The ESPResSo project
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

## This is a sample script that shows how the correlation
## module is supposed to work. It should be extended soon, but
## should already give an overview of the correlation engine.

## First set up particles in a simulation box
## with (not so important MD parameters)
set box_l 10. 
setmd box_l $box_l $box_l $box_l
set friction 1.0
set time_step 0.01
set temperature 1.0
set run_time 10000;
set int_steps 100; # number of steps per integration round

thermostat langevin $temperature $friction
setmd time_step $time_step
setmd skin 0.1
t_random seed [ pid ]

# set up some non-interacting particles of type 0
part 0 pos 0. 0. 0. type 0
part 1 pos 1. 1. 1. type 0

# Define some observables
################################

#velocities of particles of type 0
set vel [observable new particle_velocities type [list 0]]

# positions of all particles 
set pos [observable new particle_positions all]

# center of mass of particles with ids 0 and 2
set com_pos [observable new com_position id [list 0 1]]

# Particle specifications always refer to currently existing particles
# and are internally translated to a list of particle ids
# if we add more particles later, they will not be accounted for

# velocity autocorrelation function of particles of type 0
set vacf1 [correlation new obs1 $vel corr_operation scalar_product tau_max 1 dt $time_step]
# this is the minimum number of arguments to the correlation
# by default it uses the trivial correlation algorithm, which is usable only for relatively
# short tau_max < 100*dt

# same vacf as above, but using the multiple tau correlator and much longer tau_max
set vacf2 [correlation new obs1 $vel corr_operation scalar_product tau_lin 16 tau_max 100 dt $time_step compress1 linear]
# for longer time scales, use the multiple tau correlator by specifying tau_lin
# linear compression avoids loss of statistical quality of the data on long
# time scales, but is only usable with some combinations of observables 
# and correlation operations

# mean square displacement of all particles
set msd [correlation new obs1 $pos corr_operation square_distance_componentwise tau_lin 16 tau_max $run_time dt $time_step compress1 discard1]
# we want to compute the msd for for all time scales between time step and total simulation time
# therefore we use the multiple tau correlator with tau_lin=16 
# discard1 is the default compression function of the multiple tau correlator 
# to compute msd, discard1 is the only safe choice

# same msd as above, but we will update it manually, with much lower frequency
set msd_man [correlation new obs1 $pos corr_operation square_distance_componentwise tau_lin 16 tau_max $run_time dt [expr $time_step*$int_steps] compress1 discard1]

# FIXME tcl_input used to be implemented, but somehow got lost?
# msd of particle 1 based on TCL input, but with much lower sampling frequency
#set dim_input 3; # dimensionality of the input needs to be specified
#set msd_tcl [correlation new obs1 tclinput $dim_input corr_operation scalar_product tau_lin 20 tau_max 100 delta_t [expr $time_step*$int_steps] compress1 discard1]
# this will be also updated manually 

# Tell Espresso to update the desired correlations automatically
correlation $vacf1 autoupdate start
correlation $vacf2 autoupdate start
correlation $msd autoupdate start


## Now we want to measure the mobility of particle 1. We use an external force 
## and investigate the mean of its velocity in x direction.

# create one more particle of type 1 with an external force applied 
set force 1.
set new_id [setmd n_part];
part $new_id pos 1 1 1 ext_force $force 0. 0. type 1

#velocities of particles of type 0
set vel_force [observable new particle_velocities type [list 1]]

# now we take the componentwise product instead of scalar product to get 
# vacf in x, y and z separately 
set vacf_force [correlation new obs1 $vel_force corr_operation componentwise_product tau_lin 16 tau_max 100. dt $time_step compress1 linear]
# and make it autoupdate
correlation $vacf_force autoupdate start


# FIXME the rest has to be updated by Stefan, because I do not completely understand, what kind of results is should produce

#correlation [ correlation n_corr ] first_obs density_profile type { 0 } startz 0 stopz 10 nbins 10 second_obs density_profile id { 0 } startz 0 stopz 10 nbins 10 corr_operation componentwise_product tau_lin 10 tau_max 1. delta_t [ setmd time_step ] compress1 discard1 compress2 discard2
#correlation 1 first_obs radial_density_profile type { 0 } center 5. 5. 5. stopr 5 nbins 10 second_obs radial_density_profile type { 0 } center 5. 5. 5. stopr 5 nbins 10 corr_operation componentwise_product tau_lin 5 hierarchy_depth 1 delta_t [ setmd time_step ] compress1 discard1 compress2 discard2

#inter 0 0 lennard-jones 0 0.25 0. 0. 
## We also calculate the variance of the x component of the velocity 
## as reference value (to see that everything works).
set var 0.
set av 0.


## Now comes the main integration loop
set round 0;
set time [setmd time];
#set ofile [ open "v.dat" "w"]
while { $time < $run_time } {
  if { [expr $round%1000] == 0 } { 
    puts "Integration round $round, time $time";
  }
  if { [integrate $int_steps] != "" } { 
    puts "integration failed";
    exit;
  } else { incr round; }
  # Explicit call to update a correlation 
  correlation $msd_man update;
  # Updating the correlation from TCL input
  #correlation $msd_tcl update [ part 0 print v ];
  set av [ expr $av + [ lindex [ part 0 print v ] 0 ] ]
  set var [expr $var + [ lindex [ part 0 print v ] 0 ] *  [ lindex [ part 0 print v ] 0 ] ]
  #puts $ofile [ part 0 print v ]
  set time [setmd time];

}
#close $ofile
#set file_corr_number [ correlation n_corr ]
#correlation $file_corr_number first_obs textfile "v.dat" corr_operation scalar_product tau_lin 20 tau_max 100. delta_t .01 compress1 discard1
#correlation $file_corr_number update_from_file
#correlation $file_corr_number write_to_file "test.dat"

# to make use of all the history, finalize all correlations when the integration is done
for {set i 0} {$i < [correlation n_corr] } {incr i} {
  correlation $i finalize;
}

correlation $vacf1 write_to_file "vacf1.dat"
correlation $vacf2 write_to_file "vacf2.dat"
correlation $msd write_to_file "msd.dat"
correlation $msd_man write_to_file "msd_man.dat"
#correlation $msd_tcl write_to_file "msd_tcl.dat"
correlation $vacf_force write_to_file "vacf_force.dat"
exit;

#Lets look at the average velocities of the particles
#with external force
set average [ correlation 4 print average1 ]
set variance [ correlation 4 print variance1 ]
set corrtime [ correlation 4 print correlation_time ]
set stdev_mean [ correlation 4 print average_errorbars]
set true_value [ list ]
set true_correlation_time [ list ]
for { set i 0 } { $i < $part_with_force } { incr i } {
  lappend true_value [ expr $force/$friction ] 0. 0.
  lappend true_correlation_time [ expr 1/$friction ]  [ expr 1/$friction ] [ expr 1/$friction ]
}
#for { set i 0 } { $i < [ llength $average ] } { incr i } {
#  if { [ lindex $corrtime $i  ] > 0 } {
#    lappend stdev_mean [ expr sqrt( [ lindex $variance $i ] / $nsteps / [ setmd time_step ] * [ lindex $corrtime $i ] ) ]
#  } else {
#    lappend stdev_mean 0.
#  }
#}

puts "average    [ lrange $average    0 5 ] "
puts "stdev_mean [ lrange $stdev_mean 0 5 ] "
puts "true value [ lrange $true_value 0 5 ] "
puts "corrtime   [ lrange $corrtime   0 5 ] "
puts "true ctime [ lrange $true_correlation_time 0 5 ]"
set counter 0.
for { set i 0 } { $i < [ llength $average ] } { incr i } {
  if { abs([ lindex $average $i ] - [ lindex $true_value $i ]) < [ lindex $stdev_mean $i ] } {
    set counter [ expr $counter + 1 ]
  }
}
puts "[ expr $counter/([llength $average ]) *100] % are within 1 sigma"

exit

## Finally we print the result to the screen.

#puts "average [ expr $av/$nsteps ] [ lindex [ correlation 0 print average1 ] 0 ]"
#puts "variance [ expr $var/$nsteps ] [ lindex [ correlation 0 print variance1 ] 0 ]"
#set ct [ correlation 0 print correlation_time ]
#correlation 0 print
#set dens [ correlation 0 print average1 ] 
puts [ correlation 0 print average1 ] 
#puts [ correlation 1 print average1 ] 

set mean 0
set var 0
foreach { x } $ct {
  puts $x
  set mean [ expr $mean + $x ]
  set var [ expr $var + $x*$x ]
}
set mean [expr $mean / [ llength $ct ] ]
set stddev [ expr sqrt( $var/[ llength $ct ] - $mean*$mean) ]
puts "mean: $mean stddev: $stddev stderr: [ expr $stddev/sqrt([llength $ct]-1) ] "
