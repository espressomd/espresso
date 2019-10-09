from __future__ import print_function
import espressomd
from espressomd import constraints
import numpy as np
import matplotlib.pyplot as plt

espressomd.assert_features("STOKESIAN_DYNAMICS")

system = espressomd.System(box_l=[10, 10, 10])
system.time_step = 1.0
system.cell_system.skin = 0.4

system.thermostat.set_sd(viscosity=1.0, device="cpu", radii={0: 1.0})
system.integrator.set_sd()

system.part.add(pos=[-5, 0, 0], rotation=[1, 1, 1])
system.part.add(pos=[0, 0, 0], rotation=[1, 1, 1])
system.part.add(pos=[7, 0, 0], rotation=[1, 1, 1])

gravity = constraints.Gravity(g=[0, -1, 0])
system.constraints.add(gravity)

intsteps = 13000
pos = np.empty([intsteps, 3 * len(system.part)])
for i in range(intsteps):
    system.integrator.run(1)
    for n, p in enumerate(system.part):
        pos[i, 3 * n:3 * n + 3] = p.pos

for n, p in enumerate(system.part):
    plt.plot(pos[:, 3 * n], pos[:, 3 * n + 1])
plt.show()
