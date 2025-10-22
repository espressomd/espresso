#
# Copyright (C) 2013-2022 The ESPResSo project
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
import matplotlib.pyplot as plt
import numpy as np
import espressomd
import tqdm

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
        system.part.add(type=0, pos=np.random.random(3) * system.box_l)
    else:
        system.part.add(type=1, pos=np.random.random(3) * system.box_l)

print(
    f"Simulate {n_part} particles in a cubic box {box_l} at density {density}.")
print("Interactions:\n")
act_min_dist = system.analysis.min_dist()
print(f"Start with minimal distance {act_min_dist}")


#############################################################
#  Warmup Integration                                       #
#############################################################

print(f"""\
Start warmup integration:
At maximum {warm_n_times} times {warm_steps} steps
Stop if minimal distance is larger than {min_dist}""")
print(system.non_bonded_inter[0, 0].lennard_jones)

# minimize energy using min_dist as the convergence criterion
system.integrator.set_steepest_descent(f_max=0, gamma=1e-3,
                                       max_displacement=lj_sig / 100)
i = 0
while i < warm_n_times and system.analysis.min_dist() < min_dist:
    print(f"minimization: {system.analysis.energy()['total']:+.2e}")
    system.integrator.run(warm_steps)
    i += 1

print(f"minimization: {system.analysis.energy()['total']:+.2e}")
print()
system.integrator.set_vv()

# activate thermostat
system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)

# Just to see what else we may get from the C++ core
import pprint
pprint.pprint(system.cell_system.get_state(), width=1)


#############################################################
#      Integration                                          #
#############################################################
print(f"\nStart integration: run {int_n_times} times {int_steps} steps")

types = [0, 1]
sf_Sk_avg = None
for i in tqdm.trange(int_n_times):
    system.integrator.run(int_steps)
    sf_k, sf_Sk = system.analysis.structure_factor(sf_types=types, sf_order=20)
    if sf_Sk_avg is None:
        sf_Sk_avg = np.zeros(sf_Sk.shape, dtype=float)
    sf_Sk_avg += np.copy(sf_Sk) / int_n_times

plt.plot(sf_k, sf_Sk_avg)
plt.xlabel("Wavevector $q$")
plt.ylabel("Structure factor $S(q)$")
plt.show()
