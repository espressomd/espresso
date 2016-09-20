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
from espressomd import code_info

import cPickle as pickle
import os 
import numpy as np

print("""
=======================================================
=                    lj_tutorial.py                   =
=======================================================

Program Information:""")
print(code_info.features())

# System parameters
#############################################################
n_part  = 108
density = 0.8442

skin        = 0.1
time_step   = 0.001 
eq_tstep    = 0.0001
temperature = 0.728

box_l       = np.power(n_part/density, 1.0/3.0) + 2*skin

warm_steps  = 100
warm_n_time = 2000
min_dist    = 0.87

# integration
sampling_interval       = 10
equilibration_interval  = 1000

sampling_iterations     = 10000
equilibration_iterations= 20


# Interaction parameters (Lennard Jones)
#############################################################

lj_eps = 1.0
lj_sig = 1.0
lj_cut = 2.5*lj_sig
lj_cap = 20


# System setup
#############################################################
system              = espressomd.System()

if not os.path.exists('data') :
    os.mkdir('data')

system.time_step    = time_step
system.cell_system.skin         = skin

system.box_l = [box_l, box_l, box_l]

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig,
    cutoff=lj_cut, shift="auto")
system.non_bonded_inter.set_force_cap(lj_cap)

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
    system.non_bonded_inter.set_force_cap(lj_cap)

system.non_bonded_inter.set_force_cap(0)

print("\nWarm up finished\n")

system.time_step = eq_tstep 

for i in range(equilibration_iterations):
    system.integrator.run(equilibration_interval)
    energies = system.analysis.energy()
    print("eq run {} at time {}\n".format(i, system.time))

print("\nEquilibration done\n")
# Switch thermostat off NVE msd measurement
system.thermostat.turn_off()

print("\nSampling\n")
system.time_step = time_step

en_fp   = open('data/energy.dat', 'w')
sys_fp  = open('data/sim_info.pickle', 'w')
part_fp = open('data/part.pickle', 'w')
msd_fp  = open('data/msd.dat', 'w')
rdf_fp  = open('data/rdf.dat', 'w')
en_fp.write("#\n#\n#\n# Pressure   Kinetic Potential   Temperature\n#\n")

#start positions of all particles for MSD calculation on the fly
start_pos=system.part[:].pos

# save system setup
pickle.dump(system, sys_fp, -1)

msd = np.zeros((sampling_iterations,))

# Data arrays for simple error estimation
etotal = np.zeros((sampling_iterations,))
ptotal = np.zeros((sampling_iterations,))

# save start particle configuration
pickle.dump(system.part, part_fp, -1)

# analyzing the radial distribution function
# setting the parameters for the rdf
r_bins = 30
r_min  = 0.0
r_max  = system.box_l[0]/2.0

avg_rdf=np.zeros((r_bins,))

for i in range(1, sampling_iterations + 1):
    system.integrator.run(sampling_interval)
    energies = system.analysis.energy()
    pressure = system.analysis.pressure()

    r, rdf = system.analysis.rdf(rdf_type="rdf", type_list_a=[0], type_list_b=[0], r_min=r_min, r_max=r_max, r_bins=r_bins)
    avg_rdf+= rdf/sampling_iterations

    kinetic_temperature = energies['ideal']/( 1.5 * n_part)

    en_fp.write("%i\t%1.5e\t%1.5e\t%1.5e\t%1.5e\t%1.5e\n" % (i, pressure['total'], energies['total'], energies['ideal'], energies['total'] - energies['ideal'], kinetic_temperature))

    etotal[i-1] = energies['total']
    ptotal[i-1] = pressure['total']

    ####################################################
    # TODO:
    # For an efficient calculation of the MSD and VACF the correlator and
    # observable features are needed in python
    ####################################################
    for j in range(n_part):
        msd[i-1] += 1.0/n_part * ( (start_pos[j,0] - system.part[j].pos[0])**2 +
                (start_pos[j,1] - system.part[j].pos[1])**2 +
                (start_pos[j,2] - system.part[j].pos[2])**2 )
    msd_fp.write("%i %1.5e\n" % (i, msd[i-1]) )

print("\nMain sampling done\n")

# calculate the variance of the total energy and total pressure using scipys statistic operations
error_total_energy=np.sqrt(etotal.var())/np.sqrt(sampling_iterations)
error_total_pressure=np.sqrt(ptotal.var())/np.sqrt(sampling_iterations)

en_fp.write("#mean_energy energy_error mean_pressure pressure_error\n#%1.5e %1.5e %1.5e %1.5e" \
        % (etotal.mean(), error_total_energy, ptotal.mean(), error_total_pressure) )

# write out the radial distribution data
for i in range(r_bins):
    rdf_fp.write("%1.5e %1.5e\n" % (r[i], avg_rdf[i]))


sys_fp.close()
msd_fp.close()
en_fp.close()
part_fp.close()
rdf_fp.close()
