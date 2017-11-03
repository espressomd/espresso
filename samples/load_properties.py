#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
import espressomd
from espressomd import electrostatics
from espressomd import electrostatic_extensions
import numpy
try:
    import cPickle as pickle
except ImportError:
    import pickle


print("""
=======================================================
=                load_properties.py                   =
=======================================================

Program Information:""")
print(espressomd.features())

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

# Import system properties
#############################################################
system = espressomd.System()
with open("system_save", "r") as system_save:
    pickle.load(system_save)

with open("nonBondedInter_save", "r") as bond_save:
    pickle.load(bond_save)

print("Non-bonded interactions from checkpoint:")
print(system.non_bonded_inter[0, 0].lennard_jones.get_params())

print("Force cap from checkpoint:")
print(system.force_cap)


# Integration parameters
#############################################################
with open("thermostat_save", "r") as thermostat_save:
    pickle.load(thermostat_save)


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
#############################################################
# Interaction setup
#############################################################

if not system.non_bonded_inter[0, 0].lennard_jones.is_active():
    system.non_bonded_inter[0, 0].lennard_jones.set_params(
        epsilon=lj_eps, sigma=lj_sig,
        cutoff=lj_cut, shift="auto")
    system.force_cap = lj_cap
    print("Reset Lennard-Jones Interactions to:")
    print(system.non_bonded_inter[0, 0].lennard_jones.get_params())

exit()

# Import of particle properties and P3M parameters
#############################################################
with open("particle_save", "r") as particle_save:
    pickle.load(particle_save)

act_min_dist = system.analysis.mindist()

with open("p3m_save", "r") as p3m_save:
    p3m = pickle.load(p3m_save)
print(p3m.get_params())

system.actors.add(p3m)

# Check import
#############################################################
import pprint
pprint.pprint(system.cell_system.get_state(), width=1)
pprint.pprint(system.thermostat.get_state(), width=1)
# pprint.pprint(system.part.__getstate__(), width=1)
pprint.pprint(system.__getstate__(), width=1)


print("P3M parameters:\n")
p3m_params = p3m.get_params()
for key in list(p3m_params.keys()):
    print("{} = {}".format(key, p3m_params[key]))

print(system.actors)

#############################################################
#      Integration                                          #
#############################################################
print("\nStart integration: run %d times %d steps" % (int_n_times, int_steps))

# remove force capping
lj_cap = 0
system.force_cap = lj_cap
print(system.non_bonded_inter[0, 0].lennard_jones)

# print(initial energies)
energies = system.analysis.energy()
print(energies)

j = 0
for i in range(0, int_n_times):
    print("run %d at time=%f " % (i, system.time))

    system.integrator.run(int_steps)

    energies = system.analysis.energy()
    print(energies)

# terminate program
print("\nFinished.")
