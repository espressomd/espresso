from espressomd import System, lb
import numpy as np

# System setup
system = System()
system.time_step = 0.1
system.cell_system.skin = 0.2
system.box_l = [16, 16, 16]

lb_friction = 10

lbf = lb.LBFluid_GPU(agrid=1, dens=1, visc=1, tau=0.01, fric=lb_friction)
system.actors.add(lbf)
system.thermostat.set_lb(kT=1)

system.part.add(pos=[0, 0, 0])

# TODO
# # define a position observable
# set pos [observable new particle_positions all]
# # mean square displacement of all particles
# #set msd [correlation new obs1 $pos corr_operation square_distance_componentwise tau_lin 16 tau_max [expr $loops*$steps_per_loop*[setmd time_step]] dt [setmd time_step] compress1 discard1]
# set msd [correlation new obs1 $pos corr_operation square_distance_componentwise tau_lin 16 tau_max 1000.0 dt [setmd time_step] compress1 discard1]

## perform a couple of steps to come to equilbrium
print("Equilibrating the system.")
system.integrator.run(1e3)
print("Equlibration finished.")

# TODO
# correlation $msd autoupdate start

print("Sampling started.")
for i in range(4e5):
    system.integrator.run(10)

    if i % 1e3 == 0:
        sys.stdout.write("\rSampling: %05i"%i)
        sys.stdout.flush()

print("Sampling finished.")

# TODO
# correlation $msd finalize
# set msd_out [open "./msd_gamma${lb_friction}.dat" "w"]
# set number_of_datapoints [llength [correlation $msd print]]
# for { set i 0 } { $i < $number_of_datapoints } { incr i } {
# 	puts $msd_out "[lindex [lindex [correlation $msd print] $i] 0]\
#  [expr 1./3.*([lindex [lindex [correlation $msd print] $i] 2]\
#  +[lindex [lindex [correlation $msd print] $i] 3]\
#  +[lindex [lindex [correlation $msd	print] $i] 4])]";
#  }
#  close $msd_out

