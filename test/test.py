
from __future__ import print_function
import espressomd
from espressomd import thermostat
from espressomd import interactions
from espressomd import shapes

import numpy

# System parameters
#############################################################

system = espressomd.System()

#if no seed is provided espresso generates a seed

system.time_step = 0.01
system.cell_system.skin = 10.0
system.box_l = [100, 100, 100]
system.thermostat.set_langevin(kT=1.0, gamma=1.0)
system.cell_system.set_n_square(use_verlet_lists=False)

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=1, sigma=1,
    cutoff=2**(1. / 6), shift="auto")

system.part.add(id=0, pos=[1, 2, 3], type=0)

# not sure why set_parameters isn't working
#cus1=shapes.Custom()
#cus1.set_parameter("dist", -19)
#cus1.set_parameter("normal", [-1, 0, 0])

cus2=shapes.Wall(dist=-20, normal=[-1, 0, 0])
w2=system.constraints.add(particle_type=0, penetrable=0, only_positive=0, shape=cus2)



# check if the wall is actually present, the particle should be stuck near x=20

system.part[0].ext_force=[0.5, 0, 0]
for i in range(100):
    system.integrator.run(100)
#    print(system.part[0].pos)
print(system.part[0].pos)
