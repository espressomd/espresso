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


# This is a modified version of the lj liquid script, which calculates the 
# mean saure displacement versus time. 

# 1. Setup up and equilibrate the Lennard Jones liquid
# 2. Setup an observable for the particle positions
# 3. Set up a "Correlator" which records the square distance of the particles
#    for different time intervals
# 4. Integrate equations of motion
# 5. Print the result from the correlator





# 1. Setup and equilibrate LJ liquid
from __future__ import print_function
import espressomd

import os 
import numpy as np

n_part  = 200
density = 0.8442

skin        = 0.1
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


# System setup
#############################################################
system              = espressomd.System(box_l=[1.0, 1.0, 1.0])
system.seed         = system.cell_system.get_state()['n_nodes'] * [1234]
np.random.seed(system.seed)

if not os.path.exists('data') :
    os.mkdir('data')

system.time_step    = time_step
system.cell_system.skin         = skin

system.box_l = [box_l, box_l, box_l]

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig,
    cutoff=lj_cut, shift="auto")
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

print("\nWarm up finished\n")

system.time_step = eq_tstep 

for i in range(equilibration_iterations):
    system.integrator.run(equilibration_interval)
    energies = system.analysis.energy()
    print("eq run {} at time {}\n".format(i, system.time))

print("\nEquilibration done\n")

system.time_step = time_step




# 2. Setup an observable for the particle positions
from espressomd.observables import ParticlePositions
# We pass the ids of the particles to be tracked to the observable.
# Here this is 0..n_part
part_pos=ParticlePositions(ids=range(n_part))




# 3. Define a correlator recording the square distances the particles have traveled
# in various time intervals
from espressomd.correlators import Correlator
# The correlator works with the part_pos observable. Here, no second 
# observableis needed

# For short time spans, we record in linear time distances
# for larger time spans, we record in larger and larger steps.
# Use 10 linear measurements 10 time steps apart, each.

# The "square_distance_component_wise" correlation operation tells
# the correlator how to calculate a result from the measurements of the 
# observables at different times

corr=Correlator(obs1=part_pos,
                tau_lin=10,dt=10*time_step,
                tau_max=10000*time_step,
                corr_operation="square_distance_componentwise")

# We want the correlator to take measurements and calculate results 
# automatically during the integration
system.auto_update_correlators.add(corr)

# 4. Integrate the equations of motion
# Note that a longer integration would be needed to achieve good quality 
# results. This will run for ~5 minutes.
system.integrator.run(100000)


# 5. Save the result of the correlation to a file
# The 1st column contains the time, the 2nd column the number of 
# measurements, and the remaining columns the mean square displacement
# for each particle and each Cartesian component
corr.finalize()
result=corr.result()

# We average over all particls and coordinaes and store id in the 3rd column
# using numpy's averaging function
result[:,3]=np.average(result[:,3:],1)
# And save the first three columns to a file
np.savetxt("data/msd.dat",result[:,:3])

# Plot the msd (third against first column of msd.dat) with your favorite 
# plotting program and measure the diffusion coefficient by fitting
# <x^2> =2Dt

