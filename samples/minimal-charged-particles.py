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
Simulate an equal number of positively and negatively charged particles
using the P3M solver. The system is maintained at a constant temperature
using a Langevin thermostat.
"""
import espressomd

required_features = ["P3M", "WCA"]
espressomd.assert_features(required_features)

from espressomd import electrostatics
import numpy as np

# System parameters
#############################################################

box_l = 10.7437
density = 0.7

# Interaction parameters (repulsive WCA)
#############################################################

wca_eps = 1.0
wca_sig = 1.0

# Integration parameters
#############################################################
system = espressomd.System(box_l=[box_l] * 3)
np.random.seed(seed=42)

system.time_step = 0.01
system.cell_system.skin = 0.4

# warmup integration (steepest descent)
warm_steps = 20
warm_n_times = 30
# convergence criterion (particles are separated by at least 90% sigma)
min_dist = 0.9 * wca_sig

# integration
int_steps = 1000
int_n_times = 10

# Non-Bonded Interaction setup
#############################################################

system.non_bonded_inter[0, 0].wca.set_params(epsilon=wca_eps, sigma=wca_sig)

# Particle setup
#############################################################

volume = box_l**3
n_part = int(volume * density)

for i in range(n_part):
    system.part.add(id=i, pos=np.random.random(3) * system.box_l)

# Assign charge to particles
for i in range(n_part // 2 - 1):
    system.part[2 * i].q = -1.0
    system.part[2 * i + 1].q = 1.0


# Warmup
#############################################################

# minimize energy using min_dist as the convergence criterion
system.integrator.set_steepest_descent(f_max=0, gamma=1e-3,
                                       max_displacement=wca_sig / 100)
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

# P3M setup after charge assigned
#############################################################
p3m = electrostatics.P3M(prefactor=1.0, accuracy=1e-2)
system.actors.add(p3m)

#############################################################
#      Integration                                          #
#############################################################

for i in range(int_n_times):
    system.integrator.run(int_steps)

    energies = system.analysis.energy()
    print(energies)
