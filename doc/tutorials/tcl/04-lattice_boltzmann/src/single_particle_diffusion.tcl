## check if a command line argument is given
if { $argc != 0 } {
  puts "single_particle_diffusion.tcl"
  puts "-----------"
  puts "Diffusion of a single polymer chain in a lattice boltzmann fluid"
  puts "usage:"
  puts "./Espresso single_particle_diffusion.tcl"
  puts ""
  puts "output: msd_gammaX.dat"
  puts "This file contains the time dependent mean square displacement."
}

## set up the box
setmd box_l 16. 16. 16.
## the skin is not important for a single particles
setmd skin 0.2
## the MD step in this simulation could be bigger, but
## we keep it consistent with typical LJ simulations
setmd time_step 0.01
## number of integration loops
set loops 400000
## number of integration steps per loop
set steps_per_loop 10


## this creates the LB fluid with the parameters
## of our choice.
set lb_friction [expr [lindex $argv 0]]
lbfluid gpu grid 1. dens 1. visc 1. tau 0.01 friction $lb_friction
## this activates the thermalization of the
## LB fluid and adds the random forces on the 
## particle
thermostat lb 1.

## create a particle
part 0 pos 0 0 0


# define a position observable
set pos [observable new particle_positions all]
# mean square displacement of all particles
#set msd [correlation new obs1 $pos corr_operation square_distance_componentwise tau_lin 16 tau_max [expr $loops*$steps_per_loop*[setmd time_step]] dt [setmd time_step] compress1 discard1]
set msd [correlation new obs1 $pos corr_operation square_distance_componentwise tau_lin 16 tau_max 1000.0 dt [setmd time_step] compress1 discard1]



## perform a couple of steps to come to equilbrium
puts "Equilibrating the system."
integrate 1000
puts "Done."

## start the correlation
correlation $msd autoupdate start
## now the production run!
for { set i 0 } { $i < $loops } { incr i } {
  if { $i % 1000 == 0 } {
    ## tell the user we are still working
    puts -nonewline "cycle $i completed\r"; flush stdout
  }
  integrate $steps_per_loop
}

correlation $msd finalize
set msd_out [open "./msd_gamma${lb_friction}.dat" "w"]
set number_of_datapoints [llength [correlation $msd print]]
for { set i 0 } { $i < $number_of_datapoints } { incr i } {
	puts $msd_out "[lindex [lindex [correlation $msd print] $i] 0]\
 [expr 1./3.*([lindex [lindex [correlation $msd print] $i] 2]\
 +[lindex [lindex [correlation $msd print] $i] 3]\
 +[lindex [lindex [correlation $msd	print] $i] 4])]";
 }
 close $msd_out
