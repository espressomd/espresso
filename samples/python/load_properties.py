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
import espressomd._system as es
import espressomd
from espressomd import thermostat
from espressomd import code_info
from espressomd import analyze
from espressomd import integrate
from espressomd import electrostatics
from espressomd import electrostatic_extensions
import numpy
import cPickle as pickle


print("""
=======================================================
=                load_properties.py                   =
=======================================================

Program Information:""")
print(code_info.features())

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
pickle.load(open("system_save","r"))

# Integration parameters
#############################################################
thermostat.Thermostat().set_langevin(1.0, 1.0)

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
system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig,
    cutoff=lj_cut, shift="auto")
system.non_bonded_inter.set_force_cap(lj_cap)

print(system.non_bonded_inter[0, 0].lennard_jones.get_params())

# Import of particle properties and P3M parameters
#############################################################
pickle.load(open("particle_save","r"))

act_min_dist = analyze.mindist(es)

p3m=pickle.load(open("p3m_save","r"))
print(p3m.get_params())

system.actors.add(p3m)

# Check import
#############################################################
print("""
ro variables:
cell_grid     {0.cell_grid}
cell_size     {0.cell_size} 
local_box_l   {0.local_box_l} 
max_cut       {0.max_cut}
max_part      {0.max_part}
max_range     {0.max_range} 
max_skin      {0.max_skin}
n_nodes       {0.n_nodes}
n_part        {0.n_part}
n_part_types  {0.n_part_types}
periodicity   {0.periodicity}
transfer_rate {0.transfer_rate}
""".format(system))


print("P3M parameters:\n")
p3m_params = p3m.get_params()
for key in p3m_params.keys():
    print("{} = {}".format(key, p3m_params[key]))

print(system.actors)

#############################################################
#      Integration                                          #
#############################################################
print("\nStart integration: run %d times %d steps" % (int_n_times, int_steps))

# remove force capping
lj_cap = 0
system.non_bonded_inter.set_force_cap(lj_cap)
print(system.non_bonded_inter[0, 0].lennard_jones)

# print initial energies
energies = analyze.energy(system=system)
print(energies)

j = 0
for i in range(0, int_n_times):
    print("run %d at time=%f " % (i, system.time))

    integrate.integrate(int_steps)

    energies = analyze.energy(system=system)
    print(energies)

# terminate program
print("\nFinished.")
