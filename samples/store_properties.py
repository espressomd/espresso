#
# Copyright (C) 2013,2014 The ESPResSo project
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
from __future__ import print_function
import numpy
import espressomd
from espressomd import electrostatics
from espressomd import electrostatic_extensions
from samples_common import open

print("""
=======================================================
=                  store_properties.py                =
=======================================================

Program Information:""")
print(espressomd.features())
espressomd.assert_features(["ELECTROSTATICS"])

dev = "cpu"

# System parameters
#############################################################

# 10 000  Particles
box_l = 10.7437
density = 0.7

# Interaction parameters (repulsive Lennard Jones)
#############################################################

lj_eps = 1.0
lj_sig = 1.0
lj_cut = 1.12246
lj_cap = 20

# Integration parameters
#############################################################
system = espressomd.System()
system.time_step = 0.01
system.cell_system.skin = 0.4
system.thermostat.set_langevin(kT=1.0, gamma=1.0)


# warmup integration (with capped LJ potential)
warm_steps = 100
warm_n_times = 30
# do the warmup until the particles have at least the distance min__dist
min_dist = 0.9

# integration
int_steps = 1000
int_n_times = 10


#############################################################
#  Setup System                                             #
############################################################

# Interaction setup
#############################################################

system.box_l = [box_l, box_l, box_l]

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig,
    cutoff=lj_cut, shift="auto")
system.force_cap = lj_cap

print(system.non_bonded_inter[0, 0].lennard_jones.get_params())

# Particle setup
#############################################################

volume = box_l * box_l * box_l
n_part = int(volume * density)

for i in range(n_part):
    system.part.add(id=i, pos=numpy.random.random(3) * system.box_l)


print("Simulate {} particles in a cubic simulation box {} at density {}."
      .format(n_part, box_l, density).strip())
print("Interactions:\n")
act_min_dist = system.analysis.mindist()
print("Start with minimal distance {}".format(act_min_dist))

system.cell_system.max_num_cells = 2744


# Assingn charge to particles
for i in range(n_part / 2 - 1):
    system.part[2 * i].q = -1.0
    system.part[2 * i + 1].q = 1.0


# P3M setup after charge assigned
#############################################################
p3m = electrostatics.P3M(bjerrum_length=1.0, accuracy=1e-2)
system.actors.add(p3m)


#############################################################
#  Warmup Integration                                       #
#############################################################

print("""
Start warmup integration:
At maximum {} times {} steps
Stop if minimal distance is larger than {}
""".strip().format(warm_n_times, warm_steps, min_dist))

# set LJ cap
lj_cap = 20
system.force_cap = lj_cap
print(system.non_bonded_inter[0, 0].lennard_jones)


# Warmup Integration Loop
i = 0
while (i < warm_n_times and act_min_dist < min_dist):
    system.integrator.run(warm_steps)
    # Warmup criterion
    act_min_dist = system.analysis.mindist()
    i += 1

#   Increase LJ cap
    lj_cap = lj_cap + 10
    system.force_cap = lj_cap


# Just to see what else we may get from the c code
import pprint
pprint.pprint(system.cell_system.get_state(), width=1)
pprint.pprint(system.thermostat.get_state(), width=1)
# pprint.pprint(system.part.__getstate__(), width=1)
pprint.pprint(system.__getstate__(), width=1)


# Pickle data
###########################################################
try:
    import cPickle as pickle
except ImportError:
    import pickle

with open("particle_save", "w") as particle_save:
    pickle.dump(system.part, particle_save, -1)

with open("p3m_save", "w") as p3m_save:
    pickle.dump(p3m, p3m_save, -1)

with open("system_save", "w") as system_save:
    pickle.dump(system, system_save, -1)

with open("thermostat_save", "w") as thermostat_save:
    pickle.dump(system.thermostat, thermostat_save, -1)

with open("nonBondedInter_save", "w") as bond_save:
    pickle.dump(system.non_bonded_inter, bond_save, -1)

# terminate program
print("\nFinished.")
