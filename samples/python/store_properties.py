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
import espressomd._system as es
import espressomd
from espressomd import thermostat
from espressomd import code_info
from espressomd import analyze
from espressomd import integrate
from espressomd import electrostatics
from espressomd import electrostatic_extensions
import numpy

print("""
=======================================================
=                  store_properties.py                =
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

# Integration parameters
#############################################################
system = espressomd.System()
system.time_step = 0.01
system.skin = 0.4
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

system.box_l = [box_l, box_l, box_l]

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig,
    cutoff=lj_cut, shift="auto")
system.non_bonded_inter.set_force_cap(lj_cap)

print(system.non_bonded_inter[0, 0].lennard_jones.get_params())

# Particle setup
#############################################################

volume = box_l * box_l * box_l
n_part = int(volume * density)

for i in range(n_part):
    system.part.add(id=i, pos=numpy.random.random(3) * system.box_l)

analyze.distto(system, 0)

print("Simulate {} particles in a cubic simulation box {} at density {}."
      .format(n_part, box_l, density).strip())
print("Interactions:\n")
act_min_dist = analyze.mindist(es)
print("Start with minimal distance {}".format(act_min_dist))

system.max_num_cells = 2744


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
system.non_bonded_inter.set_force_cap(lj_cap)
print(system.non_bonded_inter[0, 0].lennard_jones)

# Warmup Integration Loop
i = 0
while (i < warm_n_times and act_min_dist < min_dist):
    integrate.integrate(warm_steps)
    # Warmup criterion
    act_min_dist = analyze.mindist(es)
    i += 1

#   Increase LJ cap
    lj_cap = lj_cap + 10
    system.non_bonded_inter.set_force_cap(lj_cap)

# Just to see what else we may get from the c code
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


# Pickle data
###########################################################
import cPickle as pickle

pickle.dump(system.part,open("particle_save","w"),-1)
pickle.dump(p3m,open("p3m_save","w"),-1)
pickle.dump(system,open("system_save","w"),-1)

# terminate program
print("\nFinished.")
