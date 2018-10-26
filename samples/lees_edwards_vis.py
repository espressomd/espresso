"""
Sample Lees-Edwards simulation
==============================

Simulate a DPD fluid with the hat potential and with the Lees-Edwards boundary
condition.

"""
from __future__ import print_function
import numpy as np
import argparse

import espressomd
from espressomd.visualization import openGLLive
from espressomd.interactions import HarmonicBond

required_features = ["LEES_EDWARDS", "HAT", "DPD"]
espressomd.assert_features(required_features)

parameters_help = """steady_shear: offset, oscillatory_shear: amplitude velocity, step: offset"""

parser = argparse.ArgumentParser(
    description=__doc__,
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    "protocol",
    choices=[
        "steady_shear",
        "oscillatory_shear",
        "step",
        "off"])
parser.add_argument("parameters", type=float, nargs="+", help=parameters_help)
args = parser.parse_args()

system = espressomd.System(box_l=[10, 10, 10], time_step=1e-2)
system.cell_system.skin = 0.0
system.cell_system.set_n_square()

# DPD parameters
n_part = 300
kT = 0.01
gamma = 1.0
r_cut = 0.5
F_max = 2.0

# Activate the thermostat
system.thermostat.set_dpd(kT=kT)
system.set_random_state_PRNG()

# Set up the DPD friction interaction
system.non_bonded_inter[0, 0].dpd.set_params(
    weight_function=0, gamma=gamma, r_cut=r_cut,
    trans_weight_function=0, trans_gamma=gamma, trans_r_cut=r_cut)

# Set up the repulsive interaction
system.non_bonded_inter[0, 0].hat.set_params(F_max=F_max,
                                             cutoff=r_cut)

pos = system.box_l * np.random.random((n_part, 3))
system.part.add(pos=pos)

if args.protocol == 'steady_shear':
    system.lees_edwards.set_params(type=args.protocol, velocity=args.parameters[0])
else:
    raise Exception('not impl')

# for i in range(1000):
#  system.integrator.run(steps=100)
#  system.part.writevtk("dpd_" + str(i) + ".vtk")

v = openGLLive(system)
v.run(1)
