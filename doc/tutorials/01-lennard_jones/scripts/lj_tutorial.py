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

import os 
import numpy as np

print("""
=======================================================
=                    lj_tutorial.py                   =
=======================================================

Program Information:""")
print(espressomd.features())

# System parameters
#############################################################
n_part  = 500
density = 0.002

skin        = 0.4
time_step   = 0.01 
eq_tstep    = 0.001
temperature = 0.728

box_l       = np.power(n_part/density, 1.0/3.0) 

warm_steps  = 100
warm_n_time = 2000
min_dist    = 0.87

# integration
sampling_interval       = 100
equilibration_interval  = 1000

sampling_iterations     = 100
equilibration_iterations= 5 


# Interaction parameters (Lennard Jones)
#############################################################

lj_eps = 1.0
lj_sig = 1.0
lj_cut = 2.5*lj_sig
lj_cap = 5 


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
act_min_dist = system.analysis.min_dist()
while i < warm_n_time and act_min_dist < min_dist :
    system.integrator.run(warm_steps)
    act_min_dist = system.analysis.min_dist()
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
print("\nSampling\n")
system.time_step = time_step
# Record energy versus time
en_fp   = open('data/energy.dat', 'w')

# Record radial distribution function
rdf_fp  = open('data/rdf.dat', 'w')

en_fp.write("#\n#\n#\n# Time\ttotal energy\tkinetic energy\tlennard jones energy\ttemperature\n")



# Data arrays for simple error estimation
etotal = np.zeros((sampling_iterations,))

# analyzing the radial distribution function
# setting the parameters for the rdf
r_bins = 50
r_min  = 0.0
r_max  = system.box_l[0]/2.0

avg_rdf=np.zeros((r_bins,))

for i in range(1, sampling_iterations + 1):
    system.integrator.run(sampling_interval)
    energies = system.analysis.energy()

    r, rdf = system.analysis.rdf(rdf_type="rdf", type_list_a=[0], type_list_b=[0], r_min=r_min, r_max=r_max, r_bins=r_bins)
    avg_rdf+= rdf/sampling_iterations

    kinetic_temperature = energies['kinetic']/( 1.5 * n_part)

    en_fp.write("%f\t%1.5e\t%1.5e\t%1.5e\t%1.5e\n" % (system.time, energies['total'], energies['kinetic'], energies['total'] - energies['kinetic'], kinetic_temperature))

    etotal[i-1] = energies['total']

print("\nMain sampling done\n")

# calculate the variance of the total energy using scipys statistic operations
error_total_energy=np.sqrt(etotal.var())/np.sqrt(sampling_iterations)

en_fp.write("#mean_energy energy_error %1.5e %1.5e\n" % (etotal.mean(), error_total_energy) )

# write out the radial distribution data
for i in range(r_bins):
    rdf_fp.write("%1.5e %1.5e\n" % (r[i], avg_rdf[i]))


en_fp.close()
rdf_fp.close()
