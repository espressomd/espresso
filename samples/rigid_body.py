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
Demonstrates the construction of a rigid object by means of the
``VIRTUAL_SITES_RELATIVE`` feature.
"""

import enum
import math

import numpy as np

import espressomd
required_features = ["VIRTUAL_SITES_RELATIVE", "MASS", "ROTATIONAL_INERTIA"]
espressomd.assert_features(required_features)
import espressomd.virtual_sites
import espressomd.rotation


system = espressomd.System(box_l=[10.0] * 3)
system.virtual_sites = espressomd.virtual_sites.VirtualSitesRelative()
system.time_step = 0.01
system.thermostat.set_langevin(kT=1.0, gamma=20.0, seed=42)


class ParticleTypes(enum.IntEnum):
    CENTER = enum.auto()
    BRANCH = enum.auto()


branch_len = 5
center = 0.5 * system.box_l

# Place six branches, pointing +/-x +/-y and +/-z.
# Note that we do not make the particles virtual at this point.
# The script uses center of mass an moment of inertia analysis routines
# to obtain the position and inertia moments of the central particle.
# Once a particle is made virtual, it will no longer contribute to
# observables involving mass. Virtual sites are not integrated via
# Newton's equation of motion and therefore do not have a meaningful mass.

for direction in np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]):
    for n in range(branch_len):
        system.part.add(pos=center + (n + 1) * direction,
                        type=ParticleTypes.BRANCH.value)
        system.part.add(pos=center - (n + 1) * direction,
                        type=ParticleTypes.BRANCH.value)


center_of_mass = system.analysis.center_of_mass(
    p_type=ParticleTypes.BRANCH.value)
print("Center of mass:", center_of_mass)

# if using multiple nodes, we need to change min_global_cut to the largest
# separation
max_inter = np.max(np.linalg.norm(system.part[:].pos - center_of_mass, axis=1))
system.min_global_cut = max_inter


principal_moments, principal_axes = espressomd.rotation.diagonalized_inertia_tensor(
    system.part[:].pos, system.part[:].mass)
# in this simple case, the cluster has principal axes aligned with the box
print("Principal moments: {}, principal axes tensor: {}".format(
    principal_moments, principal_axes))

# if we rotate the arms, we have to make sure that we set the quaternion of the
# center particle accordingly while setting the principal moments of inertia
AXIS = np.array([1., 0., 0.])
ANGLE = np.pi / 4.0


def rotate_vector(vector, axis, angle):
    return axis * np.dot(axis, vector) + math.cos(angle) * np.cross(
        np.cross(axis, vector), axis) + math.sin(angle) * np.cross(axis, vector)


for p in system.part:
    p.pos = rotate_vector(p.pos - center_of_mass, AXIS, ANGLE) + center_of_mass


principal_moments, principal_axes = espressomd.rotation.diagonalized_inertia_tensor(
    system.part[:].pos, system.part[:].mass)
# after rotating the whole object the principal axes changed
print("After rotating: principal moments: {}, principal axes tensor: {}".format(
    principal_moments, principal_axes))

# place center bead
p_center = system.part.add(
    pos=center_of_mass, mass=branch_len * 6 + 1, rinertia=principal_moments,
    rotation=[1, 1, 1], type=ParticleTypes.CENTER.value, quat=espressomd.rotation.matrix_to_quat(principal_axes))

# Relate the particles that make up the rigid body to the central particle.
# This will also mark them as `virtual = True`
for p in system.part.select(type=ParticleTypes.BRANCH.value):
    p.vs_auto_relate_to(p_center.id)

for frame in range(200):
    system.integrator.run(100)

print("Simulation finished")
