################################################################################
#                                                                              #
# Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project        #
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
#               Catalytic Reactions: Enhanced Diffusion Tutorial               #
#                                                                              #
################################################################################

# Set up the random seeds

t_random seed [expr ([clock clicks]*[pid])%2147483647]

# Read activity setting from the command line

if { [llength $argv] != 1 } {
    puts "Usage: Espresso $argv0 <passive/active = 0/1>"
    exit 1
}

set active  [lindex $argv 0]

if { $active != 0 && $active != 1 } {
    puts "Usage: Espresso $argv0 <passive/active = 0/1>"
    exit 1
}

# Set the parameters

set box_l  10
set radius  3.0
set csmall  0.1
set rate    1000

# Print input parameters

puts "Box length: $box_l"
puts "Colloid radius: $radius"
puts "Particle concentration: $csmall"
puts "Reaction rate: $rate"
puts "Active or Passive: $active\n"

# Create output directory

if { $active == 0 } {
    set outdir "./passive-system"
} else {
    set outdir "./active-system"
}

file mkdir $outdir

################################################################################

# Setup system parameters

set equi_steps 250
set equi_length 100

set prod_steps 2000
set prod_length 100

set dt 0.01
setmd time_step $dt

setmd box_l $box_l $box_l $box_l

setmd periodic 1 1 1
setmd skin 0.1
setmd min_global_cut [expr 1.1*$radius]

################################################################################

# Thermostat parameters

# Catalyzer is assumed to be larger, thus larger friction
set frict_trans_colloid 20.0
set frict_rot_colloid 20.0

# Particles are small and have smaller friction
set frict_trans_part 1.0
set frict_rot_part 1.0

# Temperature
set temp 1.0

################################################################################

# Set up the swimmer

## Exercise 1 ##
# Determine the initial position of the particle, which should be in
# the center of the box.

set x0pnt ...
set y0pnt ...
set z0pnt ...

set cent [setmd n_part]
part $cent pos $x0pnt $y0pnt $z0pnt type 0 temp $temp gamma $frict_trans_colloid gamma_rot $frict_rot_colloid

# Set up the particles

## Exercise 2 ##
# Above, we have set the concentration of the particles in the
# variable $csmall.  The concentration of both species of particles is
# equal.  Determine *how many* particles of one species there are.

# There are two species of equal concentration
set nB ...
set nA $nB

puts "Number of reactive A particles:  $nA"
puts "Number of reactive B particles:  $nB"

for {set i 0} { $i < $nA } {incr i} {

    set x [expr $box_l*[t_random]]
    set y [expr $box_l*[t_random]]
    set z [expr $box_l*[t_random]]

    # Prevent overlapping the colloid
    while { [expr ($x-$x0pnt)**2 + ($y-$y0pnt)**2 + ($z-$z0pnt)**2] < [expr 1.15*$radius**2] } {
        set x [expr $box_l*[t_random]]
        set y [expr $box_l*[t_random]]
        set z [expr $box_l*[t_random]]
    }
    part [setmd n_part] pos $x $y $z type 1 temp $temp gamma $frict_trans_part gamma_rot $frict_rot_part
}

for {set i 0} { $i < $nB } {incr i} {

    set x [expr $box_l*[t_random]]
    set y [expr $box_l*[t_random]]
    set z [expr $box_l*[t_random]]

    # Prevent overlapping the colloid
    while { [expr ($x-$x0pnt)**2 + ($y-$y0pnt)**2 + ($z-$z0pnt)**2] < [expr 1.15*$radius**2] } {
        set x [expr $box_l*[t_random]]
        set y [expr $box_l*[t_random]]
        set z [expr $box_l*[t_random]]
    }
    part [setmd n_part] pos $x $y $z type 2 temp $temp gamma $frict_trans_part gamma_rot $frict_rot_part
}

puts "box: [setmd box_l], npart: [setmd n_part]\n"

################################################################################

# Set up the WCA potential

## Exercise 3 ##
# Why are there two different cutoff lengths for the LJ interaction
# catalyzer/product and catalyzer/reactant?

set eps 5.0
set sig 1.0
set shift [expr 0.25]
set roff [expr $radius - 0.5*$sig]

# central and A particles
set cut [expr 2**(1/6.)*$sig]
inter 0 1 lennard-jones $eps $sig $cut $shift $roff

# central and B particles (larger cutoff)
set cut [expr 1.5*$sig]
inter 0 2 lennard-jones $eps $sig $cut $shift $roff

################################################################################

# Set up the reaction

set cat_range [expr $radius + 1.0*$sig]
set cat_rate [expr $rate]

## Exercise 4 ##
# We have read the acticity parameter from the command line into
# $active, where 0 means off and 1 means on.  When $active = 0 we can
# simply go on, but when $active = 1 we have to set up the reaction.
# Check the $active parameter and setup a reaction for the catalyzer
# of type 0 with the reactants of type 1 and products of type 2.  The
# reaction range is stored in $cat_range, the reaction rate in
# $cat_rate.  Use the number-conserving scheme by setting swap on.

...

################################################################################

# Perform warmup

set cap 1.0
set warm_length 100

## Exercise 5 ##
# Consult the User Guide for minimize_energy to find out the
# difference to warmup with explicit force-capping.
minimize_energy $cap $warm_length [expr 1.0/20.0] 0.05

################################################################################

# Enable the thermostat

## Exercise 6 ##
# Why do we enable the thermostat only after warmup?
thermostat langevin $temp $frict_trans_colloid

################################################################################

# Perform equilibration

# Integrate
for { set k 0 } { $k < $equi_steps } { incr k } {
    puts "Equilibration: $k of $equi_steps"
    integrate $equi_length
}

puts ""

################################################################################

for { set cnt 0 } { $cnt < 5 } { incr cnt } {
    # Set up the MSD calculation

    set tmax   [expr $prod_steps*$prod_length*$dt]

    set pos_id [observable new particle_positions id $cent]
    set msd    [correlation new obs1 $pos_id \
                corr_operation square_distance_componentwise \
                dt $dt tau_max $tmax tau_lin 16]
    correlation $msd autoupdate start

    ## Exercise 7a ##
    # Construct the auto-correlators for the AVACF, using the example
    # of the MSD.

    # Initialize the angular velocity auto-correlation function
    # (AVACF) correlator
    ...

    # Perform production

    # Integrate
    for { set k 0 } { $k < $prod_steps } { incr k } {
        puts "Production [expr $cnt + 1] of 5: $k of $prod_steps"
        integrate $prod_length
    }

    # Finalize the MSD and export
    correlation $msd autoupdate stop
    correlation $msd finalize
    correlation $msd write_to_file "$outdir/msd\_$cnt.dat"

    ## Exercise 7b ##
    # Finalize the angular velocity auto-correlation function (AVACF)
    # correlator and write the result to a file.
    ...
    correlation $avacf write_to_file "$outdir/avacf\_$cnt.dat"

}

################################################################################

exit 0
