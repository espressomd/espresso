#
# Copyright (C) 2013-2019 The ESPResSo project
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
Set up a Lennard-Jones fluid maintained at a fixed temperature by a
Langevin thermostat. The particles in the system are of two types:
type 0 and type 1. Type 0 particles interact with each other via a
repulsive WCA interaction. Type 1 particles neither interact with
themselves nor with type 0 particles. The spherically averaged
structure factor of particles of type 0 and type 1 is calculated
with :meth:`~espressomd.analyze.Analysis.structure_factor()`.
See :ref:`Structure factor`.
"""
import numpy as np
import espressomd

required_features = ["LENNARD_JONES"]
espressomd.assert_features(required_features)

print("""
=======================================================
=              lj_liquid_structurefactor.py           =
=======================================================
""")

# System parameters
#############################################################

box_l = 10.7437
density = 0.7

# Interaction parameters (repulsive Lennard-Jones)
#############################################################

lj_eps = 1.0
lj_sig = 1.0
lj_cut = 2.5 * lj_sig

# Integration parameters
#############################################################
system = espressomd.System(box_l=[box_l] * 3)
np.random.seed(seed=42)

system.time_step = 0.01
system.cell_system.skin = 0.4

# warmup integration (steepest descent)
warm_steps = 20
warm_n_times = 10
# convergence criterion (particles are separated by at least 90% sigma)
min_dist = 0.9 * lj_sig

# integration
int_steps = 1000
int_n_times = 20


#############################################################
#  Setup System                                             #
#############################################################

# structure factor file
structurefactor_type_list = [0, 1]
structurefactor_order = 20
structurefactor_bins = len(system.analysis.structure_factor(
    [0], structurefactor_order)[0])
structurefactor_k = np.zeros(structurefactor_bins)
structurefactor_Sk = np.zeros(structurefactor_bins)


# Interaction setup
#############################################################
system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")

print("LJ-parameters:")
print(system.non_bonded_inter[0, 0].lennard_jones.get_params())

# Particle setup
#############################################################

volume = box_l**3
n_part = int(volume * density)

for i in range(n_part):
    if i < n_part / 2.0:
        system.part.add(type=0, id=i, pos=np.random.random(3) * system.box_l)
    else:
        system.part.add(type=1, id=i, pos=np.random.random(3) * system.box_l)

print("Simulate {} particles in a cubic box {} at density {}."
      .format(n_part, box_l, density).strip())
print("Interactions:\n")
act_min_dist = system.analysis.min_dist()
print("Start with minimal distance {}".format(act_min_dist))


#############################################################
#  Warmup Integration                                       #
#############################################################

print("""
Start warmup integration:
At maximum {} times {} steps
Stop if minimal distance is larger than {}
""".strip().format(warm_n_times, warm_steps, min_dist))
print(system.non_bonded_inter[0, 0].lennard_jones)

# minimize energy using min_dist as the convergence criterion
system.integrator.set_steepest_descent(f_max=0, gamma=1e-3,
                                       max_displacement=lj_sig / 100)
i = 0
while i < warm_n_times and system.analysis.min_dist() < min_dist:
    print("minimization: {:+.2e}".format(system.analysis.energy()["total"]))
    system.integrator.run(warm_steps)
    i += 1

print("minimization: {:+.2e}".format(system.analysis.energy()["total"]))
print()
system.integrator.set_vv()

# activate thermostat
system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)

# Just to see what else we may get from the C++ core
import pprint
pprint.pprint(system.cell_system.get_state(), width=1)
pprint.pprint(system.__getstate__())


#############################################################
#      Integration                                          #
#############################################################
print("\nStart integration: run {} times {} steps"
      .format(int_n_times, int_steps))

for i in range(int_n_times):
    print("run {} at time={:.2f}".format(i, system.time))

    system.integrator.run(int_steps)

    structurefactor_k, structurefactor_Sk = system.analysis.structure_factor(
        structurefactor_type_list, structurefactor_order)

    energies = system.analysis.energy()
    print(energies['total'])
    linear_momentum = system.analysis.linear_momentum()
    print(linear_momentum)


# rescale structure factor values and write out data
structurefactor_Sk /= int_n_times
table = np.column_stack([structurefactor_k, structurefactor_Sk])
np.savetxt("pylj_liquid_structurefactor.tsv", table, delimiter='\t',
           fmt='%.5e', header="k,S(k)")
# reload: sf_k, sf_Sk = np.loadtxt("pylj_liquid_structurefactor.tsv").T

# terminate program
print("\nFinished.")
