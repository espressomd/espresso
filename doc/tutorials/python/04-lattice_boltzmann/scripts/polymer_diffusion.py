from espressomd import System, interactions, lb

# System setup
system = System()
system.box_l = [32, 32, 32]
system.cell_system.skin = 0.2

time_step = 0.1

# Lennard-Jones interaction
system.non_bonded_inter[0,0].lennard_jones.set_params(
    epsilon=1.0, sigma=1.0, r_cut=1.226,
    shift=0.25, cutoff=0)

# Fene Bond
fene = interactions.FeneBond(k=7, d_r_max=2)
system.bonded_inter.add(fene)

# TODO
# polymer 1 $nom 1. bond 0 mode PSAW

# TODO
# set com_pos [observable new com_position all]

# TODO
# set msd [correlation new obs1 $com_pos corr_operation square_distance_componentwise tau_lin 16 tau_max $run_time dt [setmd time_step] compress1 discard1]

print("Warming up the polymer chain.")
## For longer chains (>100) an extensive 
## warmup is neccessary ...
system.time_step 0.002
system.thermostat.set_langevin(kT=1.0, gamma=10)

for i in range(100):
    system.non_bonded_inter.set_force_cap(i)
    system.integrator.run(1000)

print("Warmup finished.")
system.non_bonded_inter.set_force_cap(0)
system.integrator.run(10000)
system.time_step = time_step
system.integrator.run(50000)

system.thermostat.turn_off()

system.part[:].v = [0,0,0]

lbf = lb.LBFluid(agrid=1, dens=1, visc=5, tau=time_step, fric=5)
system.actors.add(lbf)
system.thermostat.set_lb(kT=1)

print("Warming up the system with LB fluid.")
system.integrator.run(1000)
print("LB fluid warming finished.")

print("Sampling started.")
# TODO
# correlation $msd autoupdate start
for i in range(10000):
    system.integrator.run(100)
    sys.stdout.write("\rSampling: %05i"%i)
    sys.stdout.flush()

# TODO
# correlation $msd finalize
# # write correlation data to file
# set msd_out [open "./msd_nom${nom}.dat" "w"]
# set number_of_datapoints [llength [correlation $msd print]]
# for { set i 0 } { $i < $number_of_datapoints } { incr i } {
# 	puts $msd_out "[lindex [lindex [correlation $msd print] $i] 0]\
#  [expr 1./3.*([lindex [lindex [correlation $msd print] $i] 2]\
#  +[lindex [lindex [correlation $msd print] $i] 3]\
#  +[lindex [lindex [correlation $msd	print] $i] 4])]";
#  }
# close $msd_out
# set rh_out [open "./rh_nom${nom}.dat" "w"]
# puts $rh_out [analyze <rh> 0 1 $nom]
# close $rh_out
