################################################################################
#                                                                              #
# Copyright (C) 2010,2011,2012,2013,2014, 2015,2016 The ESPResSo project            #
#                                                                              #
# This file is part of ESPResSo.                                               #
#                                                                              #
# ESPResSo is free software: you can redistribute it and/or modify             #
# it under the terms of the GNU General Public License as published by         #
# the Free Software Foundation, either version 3 of the License, or            #
# (at your option) any later version.                                          #
#                                                                              #
# ESPResSo is distributed in the hope that it will be useful,                  #
# but WITHOUT ANY WARRANTY; without even the implied warranty of               #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                #
# GNU General Public License for more details.                                 #
#                                                                              #
# You should have received a copy of the GNU General Public License            #
# along with this program.  If not, see <http://www.gnu.org/licenses/>.        #
#                                                                              #
################################################################################
#                                                                              #
#                  Active Matter: Swimmer Flow Field Tutorial                  #
#                                                                              #
################################################################################

## Exercise 1 ##
# Create a routine to read in the hydrodynamic type
# (pusher/puller) and position at which the particle
# is initiated, set the variables 'type' and 'pos' to
# these values, respectively.

...

set type ...
set pos ...

################################################################################

set PI [expr acos(-1.0)]
set tcl_precision 8

## Exercise 2 ##
# Create an output directory that is labeled according
# to the value of the type and position, use the parameter
# 'outdir' to store this path

set outdir ...
file mkdir $outdir
after 250

# System parameters

set length 25.0
setmd box_l $length $length $length
setmd periodic 1 1 1

set prod_steps 1000
set prod_length 50
set dt 0.01
setmd time_step $dt
setmd skin 0.1
setmd min_global_cut 1.0

################################################################################

# Set the position of the particle

## Exercise 3 ##
# Determine the initial position of the particle, which
# should be in the center of the box, and shifted by 
# the value of 'pos' in the direction of the z-axis

set x0 ...
set y0 ...
set z0 ...

# Sphere size, mass, and moment of inertia, dipole force

set sph_size      0.5
set sph_mass      4.8
set Ixyz          4.8
set force         0.1 

## Exercise 4 ##
# Why is the sphere size set to 0.5 (this value is
# an approximation for the real value)? What happens when you
# change the mass and rotational inertia? Why is the value of 
# the force chosen to be low.

# Setup the particle particle

set cent [setmd n_part]
part $cent pos $x0 $y0 $z0 type 0 mass $sph_mass rinertia $Ixyz $Ixyz $Ixyz \
     swimming f_swim $force $type dipole_length [expr $sph_size + 0.5]

## Exercise 5 ##
# Why is the dipole_length chosen in this way?
# What happens if you make the length go to zero?
# Why does this happen?

################################################################################

# Setup the fluid (quiescent)

set agrid 1
set vskin 0.1
set frict 20.0
set visco 1.0
set densi 1.0
set temp 0.0

lbfluid gpu agrid $agrid dens $densi visc $visco \
            tau $dt friction $frict couple 3pt

## Exercise 6 ##
# What does 'couple 3pt' imply?
# Can the particle rotate in the flow field?

thermostat lb $temp

################################################################################

# Output the coordinates

set outfile [open "$outdir/trajectory.dat" "w"]
puts $outfile "####################################################"
puts $outfile "#        time        position       velocity       #"
puts $outfile "####################################################"
flush $outfile

# Production run

for { set k 0 } { $k <= $prod_steps } { incr k } {

  # Output quantities

  puts $outfile "[setmd time] [part 0 print pos v]"
  flush $outfile

  # Output 50 simulations

  if { [expr $k % ($prod_steps/50)] == 0 } {
    set num [expr $k/($prod_steps/50)]
    lbfluid print vtk velocity "$outdir/lb_velocity\_$num.vtk" 
    writevtk "$outdir/position\_$num.vtk" 0
  }

  # Integrate

  integrate $prod_length
}

close $outfile

## Exercise 7 ##
# Use the snapshots and paraview to visualize the final state.
# By appropriately choosing the initial position, you can ensure
# that the swimmer is in the center of the box. Explain why
# the flow lines look the way they do.

################################################################################

exit 0

################################################################################
