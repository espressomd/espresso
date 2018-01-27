
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


# This is a modified version of the lj liquid script, which simulates
# a two component lj liquid.
# By switching the lj interaction between the two components
# from atractive to purely repulsive, de-mixing can be achieved.

# 1. Setup and equilibrate LJ liquid
from __future__ import print_function
import espressomd

import os 
import numpy as np

from espressomd import visualization
from threading import Thread


n_part  = 200
density = 0.4442

skin        = 0.4
time_step   = 0.01 
eq_tstep    = 0.01
temperature = 0.728

box_l       = np.power(n_part/density, 1.0/3.0) 

warm_steps  = 100
warm_n_time = 2000
min_dist    = 0.87

# integration
sampling_interval       = 10
equilibration_interval  = 1000

sampling_iterations     = 10000
equilibration_iterations= 10


# Interaction parameters (Lennard Jones)
#############################################################

lj_eps = 1.0
lj_sig = 1.0
lj_cut = 2.5*lj_sig
lj_cap = 20

# This is the cutoff of the interaction between species 0 and 1.
# By setting it to 2**(1./6.) *lj_sig, it can be made purely repulsive
lj_cut_mixed =2.5 * lj_sig
lj_cut_mixed =2**(1./6.) * lj_sig


# System setup
#############################################################
system              = espressomd.System()

if not os.path.exists('data') :
    os.mkdir('data')

system.time_step    = time_step
system.cell_system.skin         = skin

system.box_l = [box_l, box_l, box_l]


# Here, lj interactions need to be setup for both components
# as well as for the mixed case of component 0 interacting with
# component 1.

# component 0
system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig,
    cutoff=lj_cut, shift="auto")

# component 1
system.non_bonded_inter[1, 1].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig,
    cutoff=lj_cut, shift="auto")

# mixed case
system.non_bonded_inter[0, 1].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig,
    cutoff=lj_cut_mixed, shift="auto")

system.force_cap = lj_cap

print("LJ-parameters:")
print(system.non_bonded_inter[0, 0].lennard_jones.get_params())

# Thermostat
system.thermostat.set_langevin(kT=temperature, gamma=1.0)

# Particle setup
#############################################################

volume = box_l * box_l * box_l

for i in range(n_part):
    system.part.add(id=i, pos=np.random.random(3) * system.box_l)
    # Every 2nd particle should be of component 1 
    if i%2==1: system.part[i].type=1

#############################################################
#  Warmup Integration                                       #
#############################################################

print("""
Start warmup integration:
At maximum {} times {} steps
Stop if minimal distance is larger than {}
""".strip().format(warm_n_time, warm_steps, min_dist))

i = 0
act_min_dist = system.analysis.mindist()
while i < warm_n_time and act_min_dist < min_dist :
    system.integrator.run(warm_steps)
    act_min_dist = system.analysis.mindist()
    print("run {} at time = {} (LJ cap= {} ) min dist = {}".strip().format(i, system.time, lj_cap, act_min_dist))
    i+=1
    lj_cap += 1.0
    system.force_cap = lj_cap

system.force_cap = 0


def loop():
    while True:
        system.integrator.run(100)
        visualizer.update()

visualizer = visualization.mayaviLive(system)
#Start simulation in seperate thread
t = Thread(target=loop)
t.daemon = True
t.start()

#Start blocking visualizer
visualizer.start()
