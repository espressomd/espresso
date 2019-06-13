# Copyright (C) 2010-2018 The ESPResSo project
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
Demonstrates the construction of a rigid object by means of the
VIRTUAL_SITES_RELATIVE feature.
"""

from __future__ import print_function

import espressomd
espressomd.assert_features(["VIRTUAL_SITES_RELATIVE"])
from espressomd import thermostat
from espressomd import integrate
from espressomd.virtual_sites import VirtualSitesRelative

import numpy as np

required_features = ["VIRTUAL_SITES_RELATIVE", "MASS", "ROTATIONAL_INERTIA"]
espressomd.assert_features(required_features)


box_l = 100
system = espressomd.System(box_l=[box_l, box_l, box_l])
system.set_random_state_PRNG()
system.seed = system.cell_system.get_state()['n_nodes'] * [1234]

system.time_step = 0.01
skin = 10.0
system.cell_system.skin = skin
system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)

### Particle types
type_centre = 0
type_A = 1


#############################################################
print("** Placing particles")
#############################################################


# start at p_id=1, and reserve p_id=0 for the centre bead.
p_id = 1

branch_len = 5
x = box_l * 0.5
y = box_l * 0.5
z = box_l * 0.5

# place six branches, pointing +/-x +/-y and +/-z
for direction in np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]):
    for n in range(branch_len):
        system.part.add(pos=[x, y, z] + (n + 1) * direction,
                        type=type_A, virtual=1)
        system.part.add(pos=[x, y, z] - (n + 1) * direction,
                        type=type_A, virtual=1)

system.virtual_sites = VirtualSitesRelative(have_velocity=True)

# here we calculate the center of mass position (com) and the moment of
# inertia (momI) of the object
com = system.analysis.center_of_mass(p_type=type_A)
print("center of mass is:", com)

# if using multiple nodes, we need to change min_global_cut to the largest
# deparation
if system.cell_system.get_state()['n_nodes'] > 1:
    max_dist = 0
    for p in system.part:
        dist = np.sum((p.pos - com)**2)
        if dist > max_dist:
            max_dist = dist
    max_dist = np.sqrt(max_dist)
    print("max dist is {}".format(max_dist))
    system.min_global_cut = max_dist

mat_I = system.analysis.moment_of_inertia_matrix(p_type=type_A)
#in this simple case, the cluster has principal axes aligned with the box
momI = [mat_I[0, 0], mat_I[1, 1], mat_I[2, 2]]
print("moment of intertia is", momI)

#place center bead
p_center = system.part.add(
    pos=com, mass=branch_len * 6 + 1, rinertia=momI,
    rotation=[1, 1, 1], type=type_centre)

# The virtual particles relate to the center one
for p in system.part:
    if p.virtual:
        p.vs_auto_relate_to(p_center.id)

for frame in range(200):
    system.integrator.run(100)

print("**Simulation finished")
