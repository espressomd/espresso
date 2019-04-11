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
This sample sets up a polymer and tests the available cell systems.
"""

from __future__ import print_function
import time
import numpy as np
import espressomd
espressomd.assert_features(["LENNARD_JONES"])
from espressomd import polymer
from espressomd import interactions


def profile():
    cs.skin = skin
    ti = time.time()
    system.integrator.run(n_steps)
    tf = time.time()
    print("\t with skin={} ran {:d} steps in {:f} seconds. steps/sec:{:f} "
          .format(skin, n_steps, tf - ti, n_steps * 1. / (tf - ti)))


system = espressomd.System(box_l=[100, 100, 100])
system.set_random_state_PRNG()
#system.seed = system.cell_system.get_state()['n_nodes'] * [1234]
cs = system.cell_system


# domain decomposition with verlet list: three equivalent commands
cs.set_domain_decomposition()
cs.set_domain_decomposition(True)
cs.set_domain_decomposition(use_verlet_lists=True)

system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)
system.time_step = 0.01


# Create a minimal polymer
system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=1, sigma=1,
    cutoff=2**(1. / 6), shift="auto")
fene = interactions.FeneBond(k=10, d_r_max=1.5)
system.bonded_inter.add(fene)
polymer.create_polymer(N_P=1, bond_length=0.97, MPC=100,
                       bond=fene, start_pos=[0, 0, 0])

n_steps = 1000

print("Testing without verlet lists...")
cs.set_domain_decomposition(use_verlet_lists=False)
for skin in np.arange(5, 15, 1):
    profile()

print("Testing with verlet lists...")
cs.set_domain_decomposition(use_verlet_lists=True)
for skin in np.arange(5, 15, 1):
    profile()

cs.set_n_square(True)
print("Testing with N-squared ...")
for skin in np.arange(5, 15, 1):
    profile()

print("Using automatic tuning...")
skin = cs.tune_skin(min_skin=0.5, max_skin=50., tol=0.5, int_steps=100)
profile()
