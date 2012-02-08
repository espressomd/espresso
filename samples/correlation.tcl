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
set friction 1.0
set force 1.
thermostat langevin 1. $friction
setmd time_step 0.01
setmd skin 0.1
t_random seed [ pid ]

part 0 pos 0. 0. 0.

## Now we set up different correlation all calculating the VACF of particle one
## They all do the same, but show different ways of how the correlation method can be applied.
## Note that the particles for the correlation must be already created!
analyze correlation 0 first_obs particle_velocities id { 0 } second_obs particle_velocities id  { 0 } corr_operation scalar_product tau_lin 20 tau_max 100 delta_t [ setmd time_step ] compress1 discard1 
analyze correlation 1 first_obs particle_velocities id { 0 } corr_operation scalar_product tau_lin 20 tau_max 100. delta_t [ setmd time_step ] compress1 discard1
analyze correlation 1 autoupdate start
analyze correlation 2 first_obs tclinput 3 corr_operation scalar_product tau_lin 20 tau_max 100. delta_t [ setmd time_step ] compress1 discard1

## We set up a second particle of which we measure the mean square displacement
part 1 pos 1. 1. 1. 
analyze correlation 3 first_obs particle_velocities  id { 1 } corr_operation square_distance tau_lin 20 tau_max 100 delta_t [ setmd time_step ] compress1 discard1 
analyze correlation 3 autoupdate start


## Now we want to measure the mobility of particle 1. We use an external force 
## and investigate the mean of its velocity in x direction.
set part_with_force 100
set counter 2 
for { set i 0 } { $i < $part_with_force } { incr i } {
  part $counter pos 1 1 1 ext_force $force 0. 0. type 1
  incr counter
}
analyze correlation 4 first_obs particle_velocities type { 1 } corr_operation componentwise_product tau_lin 20 tau_max 100. delta_t [ setmd time_step ] compress1 linear
analyze correlation 4 autoupdate start


analyze correlation [ analyze correlation n_corr ] first_obs density_profile type { 0 } startz 0 stopz 10 nbins 10 second_obs density_profile id { 0 } startz 0 stopz 10 nbins 10 corr_operation componentwise_product tau_lin 10 tau_max 1. delta_t [ setmd time_step ] compress1 discard1 compress2 discard2
#analyze correlation 1 first_obs radial_density_profile type { 0 } center 5. 5. 5. stopr 5 nbins 10 second_obs radial_density_profile type { 0 } center 5. 5. 5. stopr 5 nbins 10 corr_operation componentwise_product tau_lin 5 hierarchy_depth 1 delta_t [ setmd time_step ] compress1 discard1 compress2 discard2

inter 0 0 lennard-jones 0 0.25 0. 0. 
## We also calculate the variance of the x component of the velocity 
## as reference value (to see that everything works).
set var 0.
set av 0.
## Now comes the main integration loop
set nsteps 100
set ofile [ open "v.dat" "w"]
for { set i 0 } { $i < $nsteps } { incr i } {
  integrate 1
  ## The correlation is updated after every MD step
  analyze correlation 0 update
  analyze correlation 2 update [ part 0 print v ]
  set av [ expr $av + [ lindex [ part 0 print v ] 0 ] ]
  set var [expr $var + [ lindex [ part 0 print v ] 0 ] *  [ lindex [ part 0 print v ] 0 ] ]
  puts $ofile [ part 0 print v ]

}
close $ofile
set file_corr_number [ analyze correlation n_corr ]
analyze correlation $file_corr_number first_obs textfile "v.dat" corr_operation scalar_product tau_lin 20 tau_max 100. delta_t .01 compress1 discard1
analyze correlation $file_corr_number update_from_file
analyze correlation $file_corr_number write_to_file "test.dat"

analyze correlation 0 write_to_file "corr0.dat"
analyze correlation 1 write_to_file "corr1.dat"
analyze correlation 2 write_to_file "corr2.dat"
analyze correlation $file_corr_number write_to_file "corr3.dat"
analyze correlation 3 write_to_file "msd.dat"

#analyze correlation 3 finalize
analyze correlation 4 write_to_file "corr_with_force.dat"

#Lets look at the average velocities of the particles
#with external force
set average [ analyze correlation 4 print average1 ]
set variance [ analyze correlation 4 print variance1 ]
set corrtime [ analyze correlation 4 print correlation_time ]
set stdev_mean [ analyze correlation 4 print average_errorbars]
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

#puts "average [ expr $av/$nsteps ] [ lindex [ analyze correlation 0 print average1 ] 0 ]"
#puts "variance [ expr $var/$nsteps ] [ lindex [ analyze correlation 0 print variance1 ] 0 ]"
#set ct [ analyze correlation 0 print correlation_time ]
#analyze correlation 0 print
#set dens [ analyze correlation 0 print average1 ] 
puts [ analyze correlation 0 print average1 ] 
#puts [ analyze correlation 1 print average1 ] 

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
