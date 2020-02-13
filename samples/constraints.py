# Copyright (C) 2010-2019 The ESPResSo project
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
"""
Confine a polymer between two slabs and check that it cannot escape
them during the entire simulation.
"""

import espressomd
from espressomd import interactions
from espressomd import shapes
from espressomd import polymer

required_features = ["WCA"]
espressomd.assert_features(required_features)

import numpy as np

# System parameters
#############################################################

box_l = 50.0
system = espressomd.System(box_l=[box_l] * 3)
np.random.seed(seed=42)

system.time_step = 0.01
system.cell_system.skin = 10.0
system.cell_system.set_n_square(use_verlet_lists=False)

system.non_bonded_inter[0, 0].wca.set_params(epsilon=1, sigma=1)

num_part = 30
wall_offset = 0.1

# create random positions in the box sufficiently away from the walls
ran_pos = np.random.uniform(low=1 + wall_offset, high=box_l - 1 - wall_offset,
                            size=(num_part, 3))
system.part.add(id=np.arange(num_part),
                pos=ran_pos, type=np.zeros(num_part, dtype=int))

# bottom wall, normal pointing in the +z direction, laid at z=0.1
floor = shapes.Wall(normal=[0, 0, 1], dist=wall_offset)
c1 = system.constraints.add(
    particle_type=0, penetrable=False, only_positive=False, shape=floor)

# top wall, normal pointing in the -z direction, laid at z=49.9, since the
# normal direction points down, dist is -49.9
ceil = shapes.Wall(normal=[0, 0, -1], dist=-(box_l - wall_offset))
c2 = system.constraints.add(
    particle_type=0, penetrable=False, only_positive=False, shape=ceil)

# create stiff FENE bonds
fene = interactions.FeneBond(k=30, d_r_max=2)
system.bonded_inter.add(fene)
# start it next to the wall to test it!
start = np.array([1, 1, 1 + wall_offset])

# polymer.linear_polymer_positions will avoid violating the constraints

positions = polymer.linear_polymer_positions(n_polymers=1, beads_per_chain=50,
                                             bond_length=1.0, seed=1234,
                                             min_distance=0.9,
                                             respect_constraints=True)
for i, pos in enumerate(positions[0]):
    id = len(system.part)
    system.part.add(id=id, pos=pos)
    if i > 0:
        system.part[id].add_bond((fene, id - 1))

# Warmup
#############################################################

minimize_steps = 20
minimize_n_times = 10
min_dist = 0.9

# minimize energy using min_dist as the convergence criterion
system.integrator.set_steepest_descent(f_max=0, gamma=1e-3,
                                       max_displacement=0.01)
act_min_dist = system.analysis.min_dist()
i = 0
while (system.analysis.min_dist() < min_dist or c1.min_dist()
       < min_dist or c2.min_dist() < min_dist) and i < minimize_n_times:
    print("minimization: {:+.2e}".format(system.analysis.energy()["total"]))
    system.integrator.run(minimize_steps)
    i += 1

print("minimization: {:+.2e}".format(system.analysis.energy()["total"]))
print()
system.integrator.set_vv()

# activate thermostat
system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)

# Integration
#############################################################

int_n_times = 300
int_steps = 1000
for t in range(int_n_times):
    system.integrator.run(int_steps)
    # print the position to see if it stays within imposed constraints
    print("({:.2f}, {:.2f}, {:.2f})".format(*system.part[0].pos))
    for i in range(num_part):
        assert wall_offset < system.part[i].pos[2] < box_l - wall_offset
