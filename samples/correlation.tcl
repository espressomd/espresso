#
# Copyright (C) 2010 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
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
## should already give an idea on how the correlations are supposed to
## work. 

## First set up 100 particle in a simulation box
## with (not so important MD parameters)
set box_l 10. 
setmd box_l $box_l $box_l $box_l
thermostat langevin 1. 1.0
setmd time_step 0.1
setmd skin 0.1
for { set i 0 } { $i < 100 } { incr i } {
  part $i pos [ expr $box_l * [ t_random ] ]  [ expr $box_l * [ t_random ] ]  [ expr $box_l * [ t_random ] ]
}

## Now we initialize the correlation number 0.
## We need to pass the first and second observable 
## we want to correlate. In this case we
## choose the particle velocities as first and second observable.
## We want to calculate their componentwise product correlation,
## i.e. we get the velocity autocorrelation separate for all 
## three different directions. We want to see 20 
## linear spaced tau values with difference $time_step
## then 10 values with time difference 2*$time_step
## and 10 values with time difference 4*time_step, therefore
## choose tau_lin 20 and a hierarchy_depth of 3.
## 
## Finally we tell the correlation that the unit of time
## is $time_step, i.e. the time difference between
## sucessive updates
analyze correlation 0 first_obs particle_velocities second_obs particle_velocities corr_operation scalar_product tau_lin 20 hierarchy_depth 3 delta_t [ setmd time_step ]

## We also calculate the variance of the x component of the velocity 
## as reference value (to see that everything works).
set var 0.
## Now comes the main integration loop
set nsteps 1000
for { set i 0 } { $i < $nsteps } { incr i } {
  integrate 1
  ## The correlation is updated after every MD step
  analyze correlation 0 update
  set var [expr $var + [ lindex [ part 0 print v ] 0 ] *  [ lindex [ part 0 print v ] 0 ] ]
}
## Finally we print the result to the screen.
puts "variance [ expr $var/$nsteps ]"
analyze correlation 0 print
