from __future__ import print_function
import numpy as np
import sys 

import espressomd
from espressomd.visualization import openGLLive
from espressomd.interactions import HarmonicBond

def set_protocol(args):
  print("This is the protocol used:", args[0])
  if len(args) == 1:
    # Protocol off
    system.lees_edwards = [str(args[0])]
  elif len(args) == 2:
    # This is the step strain or steady shear case
    system.lees_edwards = [str(args[0]), float(args[1])]
  elif len(args) == 3:
    # Oscillatory case
    system.lees_edwards = [str(args[0]), float(args[1]), float(args[2])]

args = np.empty(1)

if len(sys.argv) < 2:
  print("Please provide the protocol you want to use for testing")
  sys.exit()

if len(sys.argv) >= 2:
  for arg in sys.argv:
    args = np.append(args, arg)
args = args[2:]
print(args)

system = espressomd.System(box_l=[10,10,10], time_step=1e-2)
system.cell_system.skin = 0.0
system.cell_system.set_n_square()

# DPD parameters
n_part = 300
kT=0.01
gamma=1.0
r_cut=0.5
F_max=2.0

# Activate the thermostat
system.thermostat.set_dpd(kT=kT)
system.set_random_state_PRNG()

# Set up the DPD friction interaction
system.non_bonded_inter[0,0].dpd.set_params(
  weight_function=0, gamma=gamma, r_cut=r_cut,
  trans_weight_function=0, trans_gamma=gamma, trans_r_cut=r_cut)

# Set up the repulsive interaction
system.non_bonded_inter[0,0].hat.set_params(F_max=F_max,
                                            cutoff=r_cut)

pos = system.box_l * np.random.random((n_part, 3))
system.part.add(pos=pos)

set_protocol(args)

#for i in range(1000):
#  system.integrator.run(steps=100)
#  system.part.writevtk("dpd_" + str(i) + ".vtk")

v = openGLLive(system)
v.run(1)
