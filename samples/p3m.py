#
# Copyright (C) 2013-2019 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Simulate a Lennard-Jones liquid with charges. The P3M method is used to
calculate electrostatic interactions.
"""
import numpy as np
import espressomd

required_features = ["P3M", "WCA"]
espressomd.assert_features(required_features)

from espressomd import electrostatics
import argparse

parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group()
group.add_argument("--cpu", action="store_const", dest="mode",
                   const="cpu", help="P3M on CPU", default="cpu")
group.add_argument("--gpu", action="store_const", dest="mode",
                   const="gpu", help="P3M on GPU")
args = parser.parse_args()


print("""
=======================================================
=                      p3m.py                         =
=======================================================
""")

# System parameters
#############################################################

box_l = 10
density = 0.3

# Interaction parameters (repulsive Lennard-Jones)
#############################################################

wca_eps = 10.0
wca_sig = 1.0

# Integration parameters
#############################################################
system = espressomd.System(box_l=[box_l] * 3)
np.random.seed(seed=42)

system.time_step = 0.01
system.cell_system.skin = 0.4

# warmup integration (steepest descent)
warm_steps = 20
warm_n_times = 30
# convergence criterion (particles are separated by at least 90% sigma)
min_dist = 0.9 * wca_sig

# integration
int_steps = 1000
int_n_times = 10


#############################################################
#  Setup System                                             #
#############################################################

# Interaction setup
#############################################################

system.non_bonded_inter[0, 0].wca.set_params(epsilon=wca_eps, sigma=wca_sig)

print("LJ-parameters:")
print(system.non_bonded_inter[0, 0].wca.get_params())

# Particle setup
#############################################################

volume = box_l**3
n_part = int(volume * density)

for i in range(n_part):
    system.part.add(id=i, pos=np.random.random(3) * system.box_l)

print("Simulate {} particles in a cubic box {} at density {}."
      .format(n_part, box_l, density).strip())
print("Interactions:\n")
act_min_dist = system.analysis.min_dist()
print("Start with minimal distance {}".format(act_min_dist))


# Assign charges to particles
for i in range(n_part // 2 - 1):
    system.part[2 * i].q = -1.0
    system.part[2 * i + 1].q = 1.0

# P3M setup after charge assignment
#############################################################

print("\nSCRIPT--->Create p3m\n")
if args.mode == "gpu":
    p3m = electrostatics.P3MGPU(prefactor=2.0, accuracy=1e-2)
else:
    p3m = electrostatics.P3M(prefactor=1.0, accuracy=1e-2)

print("\nSCRIPT--->Add actor\n")
system.actors.add(p3m)

print("\nSCRIPT--->P3M parameter:\n")
p3m_params = p3m.get_params()
for key in list(p3m_params.keys()):
    print("{} = {}".format(key, p3m_params[key]))

print("\nSCRIPT--->Explicit tune call\n")
p3m.tune(accuracy=1e3)

print("\nSCRIPT--->P3M parameter:\n")
p3m_params = p3m.get_params()
for key in list(p3m_params.keys()):
    print("{} = {}".format(key, p3m_params[key]))

print(system.actors)

#############################################################
#  Warmup Integration                                       #
#############################################################

print("""
Start warmup integration:
At maximum {} times {} steps
Stop if minimal distance is larger than {}
""".strip().format(warm_n_times, warm_steps, min_dist))

# minimize energy using min_dist as the convergence criterion
system.integrator.set_steepest_descent(f_max=0, gamma=1e-3,
                                       max_displacement=wca_sig / 100)
i = 0
while i < warm_n_times and system.analysis.min_dist() < min_dist:
    print("minimization: {:+.2e}".format(system.analysis.energy()["total"]))
    system.integrator.run(warm_steps)
    i += 1

print("minimization: {:+.2e}".format(system.analysis.energy()["total"]))
print()
system.integrator.set_vv()

# activate thermostat
system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)

# Just to see what else we may get from the C++ core
import pprint
pprint.pprint(system.cell_system.get_state(), width=1)
# pprint.pprint(system.part.__getstate__(), width=1)
pprint.pprint(system.__getstate__())


#############################################################
#      Integration                                          #
#############################################################
print("\nStart integration: run {} times {} steps"
      .format(int_n_times, int_steps))

for i in range(int_n_times):
    print("run {} at time={:.2f}".format(i, system.time))

    system.integrator.run(int_steps)

    energies = system.analysis.energy()
    print(energies['total'])


# terminate program
print("\nFinished.")
