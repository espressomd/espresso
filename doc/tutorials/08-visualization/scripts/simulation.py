from __future__ import print_function

from matplotlib import pyplot
from threading import Thread

import espressomd
from espressomd import visualization
import numpy

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
system = espressomd.System(box_l=[box_l, box_l, box_l])
system.seed  = system.cell_system.get_state()['n_nodes'] * [1234]
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


system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig,
    cutoff=lj_cut, shift="auto")
system.force_cap = lj_cap

# Particle setup
#############################################################

volume = box_l * box_l * box_l
n_part = int(volume * density)

for i in range(n_part):
    system.part.add(id=i, pos=numpy.random.random(3) * system.box_l)

system.analysis.dist_to(0)
act_min_dist = system.analysis.mindist()
system.cell_system.max_num_cells = 2744

#############################################################
#  Warmup Integration                                       #
#############################################################

# set LJ cap
lj_cap = 20
system.force_cap = lj_cap

# Warmup Integration Loop
i = 0
while (i < warm_n_times and act_min_dist < min_dist):
    system.integrator.run(warm_steps)
    # Warmup criterion
    act_min_dist = system.analysis.mindist()
    i += 1

#   Increase LJ cap
    lj_cap = lj_cap + 10
    system.force_cap = lj_cap

#############################################################
#      Integration                                          #
#############################################################

# remove force capping
lj_cap = 0
system.force_cap = lj_cap

def main():
    for i in range(0, int_n_times):
        print("run %d at time=%f " % (i, system.time))
        system.integrator.run(int_steps)
main()

# terminate program
print("\nFinished.")
