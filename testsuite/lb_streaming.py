import itertools
import numpy as np

import espressomd
import espressomd.lb

system = espressomd.System(box_l=[10.0]*3)
system.cell_system.skin = 0.1
system.time_step = 0.1
AGRID = 0.5
TAU = 0.1
VISC = 1./6 * AGRID * AGRID / TAU
lbf = espressomd.lb.LBFluidGPU(agrid=AGRID, tau=TAU, visc=VISC, fric=1.0, dens=1.0)
system.actors.add(lbf)

rng = range(int(system.box_l[0]/AGRID))

for i in itertools.product(rng, rng, rng):
    lbf[i].population = np.zeros(19)

pop = np.zeros(19)
pop[1] = 1
lbf[0, 0, 0].population = pop
system.integrator.run(1)
print lbf[1, 0, 0].population[1]
