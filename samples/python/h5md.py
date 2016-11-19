"""Sample for the usage of H5MD trajectory writing
in ESPResSo.
"""


import numpy as np
import espressomd # pylint: disable=import-error
from espressomd.io.writer import h5md # pylint: disable=import-error


SYS = espressomd.System()
SYS.box_l = [10.0, 10.0, 10.0]
SYS.time_step = 0.01
SYS.thermostat.set_langevin(kT=1.0, gamma=1.0)
SYS.cell_system.skin = 0.4

for i in range(1000):
    SYS.part.add(id=i, pos=np.random.random(3) * SYS.box_l, type=23)
SYS.integrator.run(steps=10)
H5_FILE = h5md.H5md(filename="bla.h5", write_pos=True, write_vel=True,
                    write_force=True, write_type=True, write_mass=False)
H5_FILE.write()
H5_FILE.close()
