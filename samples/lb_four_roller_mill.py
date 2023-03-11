#
# Copyright (C) 2021-2023 The ESPResSo project
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
Simulate a four-roller mill via slip velocity boundary conditions.
"""

import espressomd.lb
import espressomd.shapes
import espressomd.constraints
import espressomd.observables
import espressomd.math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.ticker
import itertools
import argparse
import logging
import tqdm
import sys

espressomd.assert_features(["WALBERLA"])
logging.basicConfig(level=logging.INFO, stream=sys.stdout)

parser = argparse.ArgumentParser(epilog=__doc__)
parser.add_argument("--visualizer", action="store_true", dest="visualizer",
                    help="Run the visualizer")
args = parser.parse_args()

# set up LB system
logging.info('Setting up the lattice-Boltzmann fluid')
agrid = 0.5
grid_size = np.array([64, 64, 2])
system = espressomd.System(box_l=grid_size * agrid)
system.time_step = 0.1
if args.visualizer:
    system.time_step = 0.001
system.cell_system.skin = 0.1
lb_fluid = espressomd.lb.LBFluidWalberla(
    agrid=agrid, density=0.5, kinematic_viscosity=3.2, tau=system.time_step)
system.actors.add(lb_fluid)

# set up rollers by adding tangential slip velocities to cylinders
logging.info('Setting up the rollers')
cyl_center = agrid * (grid_size // 2 + 0.5) * [1, 1, 0]
for i, j in itertools.product(range(2), range(2)):
    cyl_offset = np.array([1 + i * 0.99 - 0.51, 1 + j * 0.99 - 0.51, 0])
    cyl = espressomd.shapes.Cylinder(
        center=agrid * (grid_size // 2 + 0.5) * cyl_offset, axis=[0, 0, 1],
        length=3 * system.box_l[2], radius=14.1 * agrid, direction=1)
    if args.visualizer:
        system.constraints.add(shape=cyl)
    lb_fluid.add_boundary_from_shape(cyl)
    surface_nodes = espressomd.lb.edge_detection(
        lb_fluid.get_shape_bitmask(cyl), system.periodicity)
    tangents = espressomd.lb.calc_cylinder_tangential_vectors(
        cyl.center, lb_fluid.agrid, 0.5, surface_nodes)
    direction = 1 if (i + j) % 2 == 0 else -1
    for node, tangent in zip(surface_nodes, tangents):
        vbb = espressomd.lb.VelocityBounceBack(0.01 * direction * tangent)
        lb_fluid[node].boundary = vbb

# the system needs to be fully symmetric
mask = np.copy(lb_fluid[:, :, :].is_boundary.astype(int))
np.testing.assert_array_equal(mask, np.flip(mask, axis=0))
np.testing.assert_array_equal(mask, np.flip(mask, axis=1))
np.testing.assert_array_equal(mask, np.flip(mask, axis=2))

if args.visualizer:
    import espressomd.visualization
    visualizer = espressomd.visualization.openGLLive(
        system,
        LB_draw_velocity_plane=True,
        LB_plane_dist=0,
        LB_plane_axis=2,
        LB_vel_scale=80,
        LB_vel_radius_scale=0.05,
        LB_plane_ngrid=24,
        LB_arrow_quality=6,
        quality_constraints=48,
        camera_position=[4, 4, 50],
        background_color=[1, 1, 1],
        velocity_arrows_type_colors=[[0, 1, 0]]
    )
    visualizer.run(1)

# equilibrate the fluid
logging.info('Integration loop')
for _ in tqdm.tqdm(range(40)):
    system.integrator.run(20)

# fetch fluid and slip velocities
boundary_mask = np.squeeze(lb_fluid[:, :, 0].is_boundary.astype(bool))
quivers_boundary = []
quivers_fluid = []
for i, j in itertools.product(range(boundary_mask.shape[0]),
                              range(boundary_mask.shape[1])):
    v_fluid = lb_fluid[i, j, 0].velocity
    if boundary_mask[i, j]:
        if np.linalg.norm(v_fluid) > 1e-10:
            quivers_boundary.append([i, j, v_fluid[0], v_fluid[1]])
    else:
        quivers_fluid.append([i, j, v_fluid[0], v_fluid[1]])

# prepare canvas
logging.info('Plotting')
plt.rcParams.update({'font.size': 16})
fig1 = plt.figure()
fig2 = plt.figure()
fig3 = plt.figure()
ax1 = fig1.add_subplot(111)
ax2 = fig2.add_subplot(111)
ax3 = fig3.add_subplot(111)

# plot fluid velocity
fluid_vel = np.mean(np.linalg.norm(
    lb_fluid[:, :, :].velocity, axis=-1), axis=-1)
mask = np.ones(fluid_vel.shape) * np.nan
mask[np.nonzero(np.squeeze(lb_fluid[:, :, 0].is_boundary))] = 0
img = ax1.imshow(fluid_vel.T, origin='lower', interpolation='bilinear')
cbar = plt.colorbar(img, ax=ax1)
cbar.set_label('Fluid velocity (MD units)', rotation=90, labelpad=10)
ax1.imshow(mask.T, origin='lower', interpolation='nearest')
ax1.set_xticks([])
ax1.set_yticks([])
ax1.set(xlabel='x-axis', ylabel='y-axis')

# plot fluid velocity between the rollers
ax2.plot(agrid * np.arange(fluid_vel.shape[1]),
         np.mean(fluid_vel[31:33, :], axis=0), label='$V(x, y=L / 2)$')
ax2.set_xticks(np.arange(0, system.box_l[1] + 1, 4.0))
ax2.set(xlabel='x-axis (MD units)', ylabel='Fluid velocity (MD units)')
ax2.legend()

# plot boundary geometry
cmap = matplotlib.colors.ListedColormap(['white', 'silver', 'silver'])
cmap_bounds = [0, 1, 2]
cmap_norm = matplotlib.colors.BoundaryNorm(cmap_bounds, cmap.N)
ax3.imshow(boundary_mask.T, origin='lower', interpolation='nearest', cmap=cmap,
           norm=cmap_norm)

# add grid lines based on minor ticks
minor_locator = matplotlib.ticker.FixedLocator(np.arange(0.5, grid_size[0], 1))
ax3.xaxis.set_minor_locator(minor_locator)
ax3.yaxis.set_minor_locator(minor_locator)
ax3.tick_params(axis='both', which='minor', length=0)
ax3.grid(which='minor', color='w', linestyle='-', linewidth=1.2, zorder=2)

# remove major ticks
ax3.set_xticks([])
ax3.set_yticks([])

# add cylinder radii
# for cyl in rollers:
#     circle = plt.Circle(
#         cyl.center[:2] / agrid - agrid, cyl.radius / agrid,
#         color='r', fill=False, zorder=3)
#     ax3.add_patch(circle)

# plot velocity field
quivers_boundary = np.array(quivers_boundary)
quivers_fluid = np.array(quivers_fluid)
ax3.quiver(quivers_boundary[:, 0], quivers_boundary[:, 1], quivers_boundary[:, 2],
           quivers_boundary[:, 3], scale=.44, width=0.002, color='black',
           zorder=4, label='slip velocity')
ax3.quiver(quivers_fluid[:, 0], quivers_fluid[:, 1], quivers_fluid[:, 2],
           quivers_fluid[:, 3], scale=.44, width=0.002, color='royalblue',
           zorder=4, label='fluid velocity')
ax3.set(xlabel='x-axis', ylabel='y-axis')
ax3.legend(framealpha=1, loc='upper right')

plt.show()
