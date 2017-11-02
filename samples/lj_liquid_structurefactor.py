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
import numpy
import espressomd
from espressomd import thermostat
from espressomd import analyze
from samples_common import open

print("""
=======================================================
=              lj_liquid_structurefactor.py           =
=======================================================

Program Information:""")
print(espressomd.features())

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
system.cell_system.skin = 0.4
#es._espressoHandle.Tcl_Eval('thermostat langevin 1.0 1.0')
system.thermostat.set_langevin(kT=1.0, gamma=1.0)

# warmup integration (with capped LJ potential)
warm_steps = 100
warm_n_times = 30
# do the warmup until the particles have at least the distance min__dist
min_dist = 0.9

# integration
int_steps = 1000
int_n_times = 20


#############################################################
#  Setup System                                             #
#############################################################

# structurefactor file
structurefactor_type_list = [0, 1]
structurefactor_order = 20

structurefactor_file = open("pylj_liquid_structurefactor.dat", "w")
structurefactor_file.write("# k\tS(k)\n")
structurefactor_bins = len(system.analysis.structure_factor(
    [0], structurefactor_order)[0])
structurefactor_k = numpy.zeros(structurefactor_bins)
structurefactor_Sk = numpy.zeros(structurefactor_bins)


# Interaction setup
#############################################################

system.box_l = [box_l, box_l, box_l]

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig,
    cutoff=lj_cut, shift="auto")
system.force_cap = lj_cap

print("LJ-parameters:")
print(system.non_bonded_inter[0, 0].lennard_jones.get_params())

# Particle setup
#############################################################

volume = box_l * box_l * box_l
n_part = int(volume * density)

for i in range(n_part):
    if i < n_part / 2.0:
        system.part.add(
            type=0, id=i, pos=numpy.random.random(3) * system.box_l)
    else:
        system.part.add(
            type=1, id=i, pos=numpy.random.random(3) * system.box_l)


system.analysis.distto(0)

print("Simulate {} particles in a cubic simulation box {} at density {}."
      .format(n_part, box_l, density).strip())
print("Interactions:\n")
act_min_dist = system.analysis.mindist()
print("Start with minimal distance {}".format(act_min_dist))

system.cell_system.max_num_cells = 2744

#############################################################
#  Warmup Integration                                       #
#############################################################

# open Observable file
obs_file = open("pylj_liquid.obs", "w")
obs_file.write("# Time\tE_tot\tE_kin\tE_pot\n")

print("""
Start warmup integration:
At maximum {} times {} steps
Stop if minimal distance is larger than {}
""".strip().format(warm_n_times, warm_steps, min_dist))

# set LJ cap
lj_cap = 20
system.force_cap = lj_cap
print(system.non_bonded_inter[0, 0].lennard_jones)

# Warmup Integration Loop
i = 0
while (i < warm_n_times and act_min_dist < min_dist):
    system.integrator.run(warm_steps)
    # Warmup criterion
    act_min_dist = system.analysis.mindist()
    i += 1


#   Increase LJ cap
    lj_cap = lj_cap + 10
    system.force_cap = lj_cap

# Just to see what else we may get from the c code
import pprint
pprint.pprint(system.cell_system.get_state(), width=1)
pprint.pprint(system.__getstate__(), width=1)


# write parameter file

set_file = open("pylj_liquid.set", "w")
set_file.write("box_l %s\ntime_step %s\nskin %s\n" %
               (box_l, system.time_step, system.cell_system.skin))


#############################################################
#      Integration                                          #
#############################################################
print("\nStart integration: run %d times %d steps" % (int_n_times, int_steps))

# remove force capping
lj_cap = 0
system.force_cap = lj_cap
print(system.non_bonded_inter[0, 0].lennard_jones)

energies = system.analysis.energy()
print(energies)

j = 0
for i in range(0, int_n_times):
    print("run %d at time=%f " % (i, system.time))

    system.integrator.run(int_steps)

    structurefactor_k, structurefactor_Sk = system.analysis.structure_factor(
        structurefactor_type_list, structurefactor_order)

    energies = system.analysis.energy()
    print(energies)
    obs_file.write('{ time %s } %s\n' % (system.time, energies))
    linear_momentum = system.analysis.analyze_linear_momentum()
    print(linear_momentum)


# rescale structure factor values and write out data
structurefactor_Sk /= int_n_times

for i in range(structurefactor_bins):
    structurefactor_file.write("{0}\t{1}\n".format(
        structurefactor_k[i], structurefactor_Sk[i]))
structurefactor_file.close()

# write end configuration
end_file = open("pylj_liquid.end", "w")
end_file.write("{ time %f } \n { box_l %f }\n" % (system.time, box_l))
end_file.write("{ particles {id pos type} }")
for i in range(n_part):
    end_file.write("%s\n" % system.part[i].pos)
    # id & type not working yet

obs_file.close()
set_file.close()
end_file.close()

# terminate program
print("\nFinished.")
