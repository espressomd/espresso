import espressomd
from espressomd.io.writer import h5md 
import numpy as np

system = espressomd.System()
system.box_l = [10.0, 10.0, 10.0]
system.time_step = 0.01
system.thermostat.set_langevin(kT=1.0, gamma=1.0)
system.cell_system.skin = 0.4

for i in range(1000):
    system.part.add(id=i, pos=np.random.random(3) * system.box_l, type=23)
system.integrator.run(steps=10)
h5 = h5md.H5md(filename="bla.h5", write_pos=True, write_vel=True,
               write_force=True, write_type=True, write_mass=False)
h5.write()
h5.close()
