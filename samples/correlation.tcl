#
# Copyright (C) 2010,2012,2013 The ESPResSo project
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

## This is a sample script that shows how the core analysis is
## supposed to work. Out showcase is the calculation of the
## diffusion coefficient of noninteracting particles, by calculating the
## MSD and by applying the Green-Kubo relation that uses the
## velocity autocorrelation function (VACF)

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
setmd skin 2.0
t_random seed [ pid ]

# set up some non-interacting particles of type 0
part 0 pos 0. 0. 0. type 0
part 1 pos 1. 1. 1. type 0

################################
# Define some observables
################################

#velocities of particles of type 0
set vel [observable new particle_velocities type [list 0]]

# positions of all particles 
set pos [observable new particle_positions all]

# center of mass of particles with ids 0 and 2
# This is not used for actual analysis, but to get
# an idea how it works
set com_pos [observable new com_position id [list 0 1]]

# Only for demonstration purpose we show an observable that 
# contains the x-velocities of both particles with implemented in TCL:
proc x_vels {} {
  set v1 [ lindex [ part 0 print v ] 0 ]
  set v2 [ lindex [ part 1 print v ] 0 ]
  return [ list $v1 $v2 ]
}
# It is twodimensional an the corresponding command is "x_vels"
set x_vels_observable [ observable new tclcommand 2 "x_vels" ]

# Particle specifications always refer to currently existing particles
# and are internally translated to a list of particle ids
# if we add more particles later, they will not be accounted for.

##################################################
# Define the corresponding correlations
##################################################

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

# Example of feeding arbitrary data to the correlator (using an input file)
# We compute msd once again, using the positions which we stored in a file

# First we need to define a procedure which enables feeding 
# the values of the tcl variable "line" to an observable tclcommand
proc my_get_line { } {
    global line;
    return $line;
}
# the observable which will be used to read positions from file
set obs_pos_file [observable new tclcommand 6 my_get_line];
# one more msd, this time we use the tclcommand observable to read positions from a file and pass them to the correlator
set msd_file [correlation new obs1 $obs_pos_file corr_operation square_distance_componentwise tau_lin 16 tau_max $run_time dt [expr $time_step*$int_steps] compress1 discard1]

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

## Now comes the main integration loop
set round 0;
set time [setmd time];

set pos_file_name "xyz.dat"
set fp [ open $pos_file_name "w"]
while { $time < $run_time } {
  if { [expr $round%1000] == 0 } { 
    puts "Integration round $round, time $time";
  }
  if { [integrate $int_steps] != "" } { 
    puts "integration failed";
    exit;
  } else { incr round; }
  # write down the velocities into the file
  puts $fp [observable $pos print];
  # Explicit call to update a correlation 
  correlation $msd_man update;
  # Updating the correlation from TCL input
  set time [setmd time];
}
close $fp;


# now read-in the positions from the file and feed them to msd_file
set fp [open $pos_file_name r];
set file_data [read $fp];
close $fp;
#  Process data file
set data [split $file_data "\n"];
# the values in variable line are used to compute observable pos_file
foreach line $data {
    if { [llength $line] == "0" } { 
        break; 
    }
    correlation $msd_file update; 
}


# to make use of all the history, finalize all correlations when the integration is done
for {set i 0} {$i < [correlation n_corr] } {incr i} {
  correlation $i finalize;
}

correlation $vacf1 write_to_file "vacf1.dat"
correlation $vacf2 write_to_file "vacf2.dat"
correlation $msd write_to_file "msd.dat"
correlation $msd_man write_to_file "msd_man.dat"
correlation $msd_file write_to_file "msd_file.dat"
correlation $vacf_force write_to_file "vacf_force.dat"
