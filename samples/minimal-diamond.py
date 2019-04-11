#
# Copyright (C) 2013-2018 The ESPResSo project
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
This sample sets up a diamond-structured polymer network.
"""
from __future__ import print_function
import espressomd
espressomd.assert_features(["LENNARD_JONES"])
from espressomd import thermostat
from espressomd import interactions
from espressomd import diamond
from espressomd.io.writer import vtf  # pylint: disable=import-error

import numpy as np

# System parameters
#############################################################

system = espressomd.System(box_l=[1.0, 1.0, 1.0])
system.set_random_state_PRNG()
system.seed = system.cell_system.get_state()['n_nodes'] * [1234]
np.seed = system.seed
system.time_step = 0.01
system.cell_system.skin = 0.4
system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)
system.cell_system.set_n_square(use_verlet_lists=False)

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=1, sigma=1,
    cutoff=2**(1. / 6), shift="auto")

system.non_bonded_inter[0, 1].lennard_jones.set_params(
    epsilon=1, sigma=1,
    cutoff=2**(1. / 6), shift="auto")

system.non_bonded_inter[1, 1].lennard_jones.set_params(
    epsilon=1, sigma=1,
    cutoff=2**(1. / 6), shift="auto")

fene = interactions.FeneBond(k=30, d_r_max=1.5)
system.bonded_inter.add(fene)

# The call to diamond.Diamond() creates 16 connected polymers.
# These polymers are initialized in a straight line connected to crosslink nodes
# Furthermore they are connected to one another across simulation boxes in a periodic fashion.
# It is crucial that the simulation box, chain length and a-parameters be
# chosen such that the arrangement will not break bonds.

# monomers per chain
MPC = 15

# length for Kremer-Grest chain
bond_length = 0.966

# The physical distance beween nodes such that a line of monomers "fit" needs to be worked out.
# This is done via the unit diamond lattice size parameter "a".
a = (MPC + 1) * bond_length / (0.25 * np.sqrt(3))

# Lastly, the created periodic connections requires a specific simulation box.
system.box_l = [a, a, a]
print("box now at ", system.box_l)

# We can now call diamond to place the monomers, crosslinks and bonds.
diamond.Diamond(a=a, bond_length=bond_length, MPC=MPC)


#############################################################
#      Warmup                                               #
#############################################################

print("Warming up...")
warm_steps = 10
lj_cap = 1
system.force_cap = lj_cap
act_min_dist = system.analysis.min_dist()

# warmp with zero temperature to remove overlaps
system.thermostat.set_langevin(kT=0.0, gamma=1.0)

# slowly ramp up the cap
while (lj_cap < 5):
    system.integrator.run(warm_steps)
    system.part[:].v = [0, 0, 0]
    lj_cap = lj_cap * 1.1
    system.force_cap = lj_cap

# remove force cap
lj_cap = 0
system.force_cap = lj_cap
system.integrator.run(warm_steps * 10)

# restore simulation temperature
system.thermostat.set_langevin(kT=1.0, gamma=1.0)
system.integrator.run(warm_steps * 10)
print("Finished warmup")


#############################################################
#      Integration                                          #
#############################################################

sim_steps = 100
print("simulating...")
system.integrator.run(sim_steps)

# because stretched polymers are not too impressive...
print("simulating a slow compression...")
for d in np.arange(1, 15):
    system.change_volume_and_rescale_particles(
        d_new=system.box_l[0] - 1, dir='xyz')
    print("box now at ", system.box_l)
    system.integrator.run(sim_steps)

# visualize at a smaller box size
outfile = open('diamond.vtf', 'w')
vtf.writevsf(system, outfile)
vtf.writevcf(system, outfile)
t_steps = 100
for t in range(t_steps):
    print("step {} of {}".format(t, t_steps))
    system.integrator.run(sim_steps)
    vtf.writevcf(system, outfile)
outfile.close()
