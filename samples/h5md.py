"""Sample for the usage of H5MD trajectory writing
in ESPResSo.
"""


import numpy as np
import espressomd  # pylint: disable=import-error
from espressomd.io.writer import h5md  # pylint: disable=import-error
from espressomd import polymer
from espressomd import interactions


sys = espressomd.System()
sys.box_l = [100.0, 100.0, 100.0]
sys.time_step = 0.01
sys.thermostat.set_langevin(kT=1.0, gamma=1.0)
sys.cell_system.skin = 0.4

fene = interactions.FeneBond(k=10, d_r_max=2)
sys.bonded_inter.add(fene)
polymer.create_polymer(N_P=5, bond_length=1.0, MPC=50, bond=fene)

sys.integrator.run(steps=0)
h5_file = h5md.H5md(filename="sample.h5", write_pos=True, write_vel=True,
                    write_force=True, write_species=True, write_mass=False, write_charge=True, write_ordered=True)
for i in range(1):
    h5_file.write()
h5_file.flush()
h5_file.close()
