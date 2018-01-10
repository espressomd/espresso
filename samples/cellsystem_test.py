#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
import espressomd as es
from espressomd import polymer
from espressomd import interactions
import time
import numpy as np


def profile():
    cs.skin = skin
    ti = time.time()
    S.integrator.run(n_steps)
    tf = time.time()
    print("\t with skin={} ran {:d} steps in {:f} seconds. steps/sec:{:f} ".format(skin,
                                                                                   n_steps, tf - ti, n_steps * 1. / (tf - ti)))


S = es.System()
cs = S.cell_system


# domain decomposition with verlet list: three equivalent commands
cs.set_domain_decomposition()
cs.set_domain_decomposition(True)
cs.set_domain_decomposition(use_verlet_lists=True)

S.thermostat.set_langevin(kT=1.0, gamma=1.0)
S.box_l = [100, 100, 100]
S.time_step = 0.01


# Create a minimal polymer
S.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=1, sigma=1,
    cutoff=2**(1. / 6), shift="auto")
fene = interactions.FeneBond(k=10, d_r_max=1.5)
S.bonded_inter.add(fene)
polymer.create_polymer(N_P=1, bond_length=0.97, MPC=100, bond=fene)

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
skin=cs.tune_skin(min_skin=0.5, max_skin=50., tol=0.5, int_steps=100)
profile()
