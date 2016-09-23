import espressomd
import sys
import numpy as np
from espressomd.io.writer.h5md import H5md

system = espressomd.System()
system.box_l = [10.0, 10.0, 10.0]
system.time_step = 0.01
system.thermostat.set_langevin(kT=1.0, gamma=1.0)
system.cell_system.skin = 0.4
for i in range(10):
    system.part.add(id=i, pos=np.random.random(3) * system.box_l, type=23)
h5md = H5md("data.h5")
system.integrator.run(steps=1)
h5md.write("p")
system.integrator.run(steps=1)
system.integrator.run(steps=1)
h5md.write("v")
system.integrator.run(steps=1)
h5md.write("pvfmt")
h5md.close()
