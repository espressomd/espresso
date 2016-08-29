from __future__ import print_function
import espressomd
from espressomd import visualization
import numpy
from matplotlib import pyplot
from threading import Thread

# System parameters
#############################################################

# 10 000  Particles
box_l = 10.7437
density = 0.7

# Interaction parameters (repulsive Lennard Jones)
#############################################################

lj_eps = 1.0
lj_sig = 1.0
lj_cut = 1.12246
lj_cap = 20

# Integration parameters
#############################################################
system = espressomd.System()
system.time_step = 0.01
system.cell_system.skin = 0.4
system.thermostat.set_langevin(kT=1.0, gamma=1.0)

# warmup integration (with capped LJ potential)
warm_steps = 100
warm_n_times = 30
# do the warmup until the particles have at least the distance min_dist
min_dist = 0.9

# integration
int_steps = 1000
int_n_times = 100


#############################################################
#  Setup System                                             #
#############################################################

# Interaction setup
#############################################################

system.box_l = [box_l, box_l, box_l]

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig,
    cutoff=lj_cut, shift="auto")
system.non_bonded_inter.set_force_cap(lj_cap)

# Particle setup
#############################################################

volume = box_l * box_l * box_l
n_part = int(volume * density)

for i in range(n_part):
    system.part.add(id=i, pos=numpy.random.random(3) * system.box_l)

system.analysis.distto(0)
act_min_dist = system.analysis.mindist()
system.cell_system.max_num_cells = 2744

#############################################################
#  Warmup Integration                                       #
#############################################################

# set LJ cap
lj_cap = 20
system.non_bonded_inter.set_force_cap(lj_cap)

# Warmup Integration Loop
i = 0
while (i < warm_n_times and act_min_dist < min_dist):
    system.integrator.run(warm_steps)
    # Warmup criterion
    act_min_dist = system.analysis.mindist()
    i += 1

#   Increase LJ cap
    lj_cap = lj_cap + 10
    system.non_bonded_inter.set_force_cap(lj_cap)

#############################################################
#      Integration                                          #
#############################################################

# remove force capping
lj_cap = 0
system.non_bonded_inter.set_force_cap(lj_cap)

energies = numpy.empty((int_steps,2))
current_time = -1
pyplot.xlabel("time")
pyplot.ylabel("energy")
plot, = pyplot.plot([0],[0])
pyplot.show(block=False)
def update_plot():
    if current_time < 0:
        return
    i = current_time
    plot.set_xdata(energies[:i+1,0])
    plot.set_ydata(energies[:i+1,1])
    pyplot.xlim(0, energies[i,0])
    pyplot.ylim(energies[:i+1,1].min(), energies[:i+1,1].max())
    pyplot.draw()
def main():
    global current_time
    for i in range(0, int_n_times):
        print("run %d at time=%f " % (i, system.time))
        system.integrator.run(int_steps)
        energies[i] = (system.time, system.analysis.energy()['total'])
        current_time = i
        update_plot()
main()

print ("Average energy: %6g" % energies[:,1].mean())

# terminate program
print("\nFinished.")
