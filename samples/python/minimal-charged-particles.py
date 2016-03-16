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
import espressomd
from espressomd import thermostat
from espressomd import analyze
from espressomd import integrate
from espressomd import electrostatics
import numpy

# System parameters
#############################################################

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
system.box_l = [box_l, box_l, box_l]
thermostat.Thermostat().set_langevin(1.0, 1.0)

# warmup integration (with capped LJ potential)
warm_steps = 100
warm_n_times = 30
# do the warmup until the particles have at least the distance min__dist
min_dist = 0.9

# integration
int_steps = 1000
int_n_times = 10


# Non-Bonded Interaction setup
#############################################################


system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig,
    cutoff=lj_cut, shift="auto")
system.non_bonded_inter.set_force_cap(lj_cap)


# Particle setup
#############################################################

volume = box_l * box_l * box_l
n_part = int(volume * density)

for i in range(n_part):
    system.part.add(id=i, pos=numpy.random.random(3) * system.box_l)

# Assingn charge to particles
for i in range(n_part / 2 - 1):
    system.part[2 * i].q = -1.0
    system.part[2 * i + 1].q = 1.0


# Warmup
#############################################################

lj_cap = 20
system.non_bonded_inter.set_force_cap(lj_cap)
i = 0
act_min_dist = analyze.mindist(system)
while (i < warm_n_times and act_min_dist < min_dist):
    integrate.integrate(warm_steps)
    # Warmup criterion
    act_min_dist = analyze.mindist(system)
    i += 1
    lj_cap = lj_cap + 10
    system.non_bonded_inter.set_force_cap(lj_cap)

lj_cap = 0
system.non_bonded_inter.set_force_cap(lj_cap)

# P3M setup after charge assigned
#############################################################
p3m = electrostatics.P3M(bjerrum_length=1.0, accuracy=1e-2)
system.actors.add(p3m)


#############################################################
#      Integration                                          #
#############################################################

for i in range(0, int_n_times):
    integrate.integrate(int_steps)

    energies = analyze.energy(system=system)
    print energies
