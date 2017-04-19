from __future__ import print_function

import espressomd
import numpy as np
from espressomd import interactions
from espressomd import polymer

system = espressomd.System()

print( espressomd.code_info.features() )

system.time_step = 0.01

# tuning parameter for frequency of Verlet rebuilds
system.cell_system.skin = 0.4
system.thermostat.set_langevin(kT=1.0, gamma=1.0)
system.cell_system.set_n_square(use_verlet_lists=False)

box_l = 100
system.box_l = [box_l, box_l, box_l]

harmonic_bond = interactions.HarmonicBond(k=100., r_0=1.0)
system.bonded_inter.add(harmonic_bond)

polymer.create_polymer(N_P = 1, bond_length = 1.0, MPC=10, start_id=0, type_poly_neutral=0, type_poly_charged=0, bond=harmonic_bond,  mode=1)
polymer.create_polymer(N_P = 1, bond_length = 1.0, MPC=10, start_id=100, type_poly_neutral=1, type_poly_charged=1, bond=harmonic_bond,  mode=1)

# add weird bond with other polymer
# this is a bond between two different particle types
system.part[5].add_bond((harmonic_bond,105))


fp_0=open('type_0.vtf', 'w')
fp_all=open('type_all.vtf', 'w')

system.part.writevsf(fp_0, types=0)
system.part.writevsf(fp_all, types='all')

for i in range(20):
    system.integrator.run(100)
    system.part.writevcf(fp_0, types=0)
    system.part.writevcf(fp_all, types='all')


exit()

