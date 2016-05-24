import sys
import espressomd._system as es
import numpy as np
from espressomd.interactions import *
from espressomd import analyze
import espressomd

#from espressomd import minimize_energy


system = espressomd.System()

box_l = 4.0 
density = 0.7

# Interaction parameters (repulsive Lennard Jones)
#############################################################

lj_eps = 1.0
lj_sig = 1.0
lj_cut = 1.12246

system.time_step = 0.01
system.skin = 0.4


system.box_l = [box_l, box_l, box_l]

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig,
    cutoff=lj_cut, shift="auto")


minimize=system.minimize_energy
minimize.init(f_max=0.0, gamma=0.01, max_steps=10000, max_displacement=0.01)

volume = box_l * box_l * box_l
n_part = int(volume * density)

for i in range(n_part):
    system.part.add(id=i, pos=np.random.random(3) * system.box_l)

minimize.minimize()

energy = analyze.energy("total") 
#print energy
if energy["total"] > 0 :
    print "energy $energy too big."
    sys.exit(1)

sys.exit(0)




