"""Sample for the usage of H5MD trajectory writing
in ESPResSo.
"""


import numpy as np
import espressomd # pylint: disable=import-error
from espressomd.io.writer import h5md # pylint: disable=import-error


sys = espressomd.System()
sys.box_l = [10.0, 10.0, 10.0]
sys.time_step = 0.01
sys.thermostat.set_langevin(kT=1.0, gamma=1.0)
sys.cell_system.skin = 0.4

for i in range(1000):
    sys.part.add(id=i, pos=np.random.random(3) * sys.box_l, type=23)
sys.integrator.run(steps=10)
h5_file = h5md.H5md(filename="bla.h5", write_pos=True, write_vel=True,
                    write_force=True, write_type=True, write_mass=False)
h5_file.write()
h5_file.close()
