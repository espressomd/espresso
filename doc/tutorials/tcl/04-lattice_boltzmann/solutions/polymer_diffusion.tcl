## check if a command line argument is given
if { $argc != 1 } {
  puts "polymer.tcl"
  puts "-----------"
  puts "Diffusion of a single polymer chain in a lattice boltzmann fluid"
  puts "Please pass the number of monomers as a command line argument!"
  puts "./Espresso polymer.tcl #monomers"
  puts ""
  puts "output: pos.dat v.dat"
  puts "The files contain the time elapsed in the simulation and the center"
  puts "of mass position and velocity"
  exit
}
set nom [ lindex $argv 0 ]
set time_step 0.01
set loops 10000
set time_steps_per_loop 100
set run_time [expr $loops*$time_steps_per_loop*$time_step]
set vmd no

setmd box_l 32. 32. 32.
setmd skin 0.2
setmd time_step 0.01

inter 0 0 lennard-jones 1. 1. 1.226 0.25 0.
inter 0 fene 7. 2.

polymer 1 $nom 1. bond 0 mode PSAW

# create an observable for the center of mass of the polymer
set com_pos [observable new com_position all]

set msd [correlation new obs1 $com_pos corr_operation square_distance_componentwise tau_lin 16 tau_max $run_time dt [setmd time_step] compress1 discard1]

if { $vmd=="yes" } {
  prepare_vmd_connection "test.psf"
  after 2000
  imd positions
} 
  
puts "Warming up the polymer chain."
## For longer chains (>100) an extensive 
## warmup is neccessary ...
setmd time_step 0.002
thermostat langevin 1. 10.
for { set i 0 } { $i < 100 } { incr i } {
  inter forcecap $i
  imd positions
  integrate 1000
}
puts "Done."
inter forcecap 0.
integrate 10000
setmd time_step $time_step
integrate 50000

thermostat off
kill_particle_motion
lbfluid gpu agrid 1. dens 1. visc 5. tau 0.01 friction 5.
thermostat lb 1.

puts "Warming up the system with LB fluid."
integrate 1000
puts "Done."
correlation $msd autoupdate start
for { set i 0 } { $i < $loops } { incr i } {
    imd positions
    integrate $time_steps_per_loop
	analyze append
	puts -nonewline "Loop $i/$loops\r"; flush stdout
}
puts "\n"
correlation $msd finalize
# write correlation data to file
set msd_out [open "./msd_nom${nom}.dat" "w"]
set number_of_datapoints [llength [correlation $msd print]]
for { set i 0 } { $i < $number_of_datapoints } { incr i } {
	puts $msd_out "[lindex [lindex [correlation $msd print] $i] 0]\
 [expr 1./3.*([lindex [lindex [correlation $msd print] $i] 2]\
 +[lindex [lindex [correlation $msd print] $i] 3]\
 +[lindex [lindex [correlation $msd	print] $i] 4])]";
 }
close $msd_out
set rh_out [open "./rh_nom${nom}.dat" "w"]
puts $rh_out [analyze <rh> 0 1 $nom]
close $rh_out
