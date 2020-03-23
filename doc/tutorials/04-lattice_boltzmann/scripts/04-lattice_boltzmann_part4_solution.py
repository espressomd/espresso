# Copyright (C) 2019 The ESPResSo project
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

import matplotlib.pyplot as plt
import numpy as np

import espressomd
espressomd.assert_features(['LB_BOUNDARIES_GPU'])
import espressomd.lb
import espressomd.lbboundaries
import espressomd.shapes

# System setup
BOX_L = 16.0
AGRID = 0.5
VISCOSITY = 2.0
FORCE_DENSITY = [0.0, 0.001, 0.0]
DENSITY = 1.5
TIME_STEP = 0.01
system = espressomd.System(box_l=[BOX_L] * 3)
system.time_step = TIME_STEP
system.cell_system.skin = 0.4

lbf = espressomd.lb.LBFluidGPU(
    agrid=AGRID, dens=DENSITY, visc=VISCOSITY, tau=TIME_STEP,
    ext_force_density=FORCE_DENSITY)
system.actors.add(lbf)

# Setup boundaries
WALL_OFFSET = AGRID
top_wall = espressomd.lbboundaries.LBBoundary(
    shape=espressomd.shapes.Wall(normal=[1, 0, 0], dist=WALL_OFFSET))
bottom_wall = espressomd.lbboundaries.LBBoundary(
    shape=espressomd.shapes.Wall(normal=[-1, 0, 0], dist=-(BOX_L - WALL_OFFSET)))

system.lbboundaries.add(top_wall)
system.lbboundaries.add(bottom_wall)

# Iterate until the flow profile converges (5000 LB updates)
system.integrator.run(5000)

# Extract fluid velocity along the x-axis
fluid_velocities = np.zeros((lbf.shape[0], 2))
for x in range(lbf.shape[0]):
    # Average over the node in y direction
    v_tmp = np.zeros(lbf.shape[1])
    for y in range(lbf.shape[1]):
        v_tmp[y] = lbf[x, y, 0].velocity[1]
    fluid_velocities[x, 0] = (x + 0.5) * AGRID
    fluid_velocities[x, 1] = np.average(v_tmp)               


def poiseuille_flow(x, force_density, dynamic_viscosity, height):
    return force_density / (2.0 * dynamic_viscosity) * \
        (height**2.0 / 4.0 - x**2.0)


# Note that the LB viscosity is not the dynamic viscosity but the
# kinematic viscosity (mu=LB_viscosity * density)
x_values = np.linspace(0.0, BOX_L, lbf.shape[0])
HEIGHT = BOX_L - 2.0 * AGRID
# analytical curve
y_values = poiseuille_flow(x_values - (HEIGHT / 2.0 + AGRID), FORCE_DENSITY[1],
                           VISCOSITY * DENSITY, HEIGHT)
# velocity is zero outside the walls
y_values[np.nonzero(x_values < WALL_OFFSET)] = 0.0
y_values[np.nonzero(x_values > BOX_L - WALL_OFFSET)] = 0.0

plt.plot(x_values, y_values, 'o-', label='analytical')
plt.plot(fluid_velocities[:, 0], fluid_velocities[:, 1], label='simulation')
plt.legend()
plt.show()
