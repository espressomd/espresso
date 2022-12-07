#
# Copyright (C) 2021 The ESPResSo project
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
Simulate a rotating cylinder in a fluid via slip velocity boundary conditions.
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

espressomd.assert_features(["WALBERLA"])

parser = argparse.ArgumentParser(epilog=__doc__)
parser.add_argument("--visualizer", action="store_true", dest="visualizer",
                    help="Run the visualizer")
args = parser.parse_args()

# set up LB system
agrid = 0.5
grid_size = np.array([31, 31, 4])
system = espressomd.System(box_l=grid_size * agrid)
system.time_step = 0.1
if args.visualizer:
    system.time_step = 0.001
system.cell_system.skin = 0.1
system.periodicity = [False, False, True]
lb_fluid = espressomd.lb.LBFluidWalberla(
    agrid=agrid, density=0.5, viscosity=3.2, tau=system.time_step)
system.actors.add(lb_fluid)

# set up cylinders
cyl_center = agrid * (grid_size // 2 + 0.5) * [1, 1, 0]
cylinder_in = espressomd.shapes.Cylinder(
    center=cyl_center, axis=[0, 0, 1], length=3 * system.box_l[2],
    radius=8.1 * agrid, direction=1)
cylinder_out = espressomd.shapes.Cylinder(
    center=cyl_center, axis=[0, 0, 1], length=3 * system.box_l[2],
    radius=14.5 * agrid, direction=-1)
lb_fluid.add_boundary_from_shape(cylinder_in)
lb_fluid.add_boundary_from_shape(cylinder_out)

# the system needs to be fully symmetric
mask = np.copy(lb_fluid[:, :, :].is_boundary.astype(int))
np.testing.assert_array_equal(mask, np.flip(mask, axis=0))
np.testing.assert_array_equal(mask, np.flip(mask, axis=1))
np.testing.assert_array_equal(mask, np.flip(mask, axis=2))

# the system needs to be closed in the x and y directions
np.testing.assert_array_equal(mask[0, :, :], 1)
np.testing.assert_array_equal(mask[-1, :, :], 1)
np.testing.assert_array_equal(mask[:, 0, :], 1)
np.testing.assert_array_equal(mask[:, -1, :], 1)

# add tangential slip velocity to the inner cylinder
velocity_magnitude = 0.01
surface_nodes = espressomd.lb.edge_detection(
    lb_fluid.get_shape_bitmask(cylinder_in), system.periodicity)
tangents = espressomd.lb.calc_cylinder_tangential_vectors(
    cylinder_in.center, lb_fluid.agrid, 0.5, surface_nodes)
for node, tangent in zip(surface_nodes, tangents):
    lb_fluid[node].boundary = espressomd.lb.VelocityBounceBack(
        velocity_magnitude * tangent)

if args.visualizer:
    import espressomd.visualization
    visualizer = espressomd.visualization.openGLLive(
        system,
        LB_draw_velocity_plane=True,
        LB_plane_dist=0,
        LB_plane_axis=2,
        LB_vel_scale=80,
        LB_vel_radius_scale=0.05,
        LB_plane_ngrid=15,
        quality_constraints=128,
        camera_position=[8, 8, 30],
        background_color=[1, 1, 1],
        velocity_arrows_type_colors=[[0, 1, 0]]
    )
    system.constraints.add(shape=cylinder_in)
    system.constraints.add(shape=cylinder_out)
    system.integrator.run(100)
    visualizer.run(1)

# add observable for the fluid velocity in cylindrical coordinates
cyl_transform_params = espressomd.math.CylindricalTransformationParameters(
    center=cyl_center, axis=[0, 0, 1], orientation=[1, 0, 0])
observable = espressomd.observables.CylindricalLBVelocityProfile(
    transform_params=cyl_transform_params,
    n_r_bins=grid_size[0] // 2,
    n_phi_bins=1,
    n_z_bins=1,
    min_r=0.0,
    max_r=system.box_l[0] / 2,
    min_phi=0.,
    max_phi=2 * np.pi,
    min_z=-system.box_l[2] / 2,
    max_z=+system.box_l[2] / 2,
    axis=[0.0, 0.0, 1.0],
    sampling_density=1
)
obs_data_baseline = observable.calculate()

# equilibrate the fluid
system.integrator.run(100)
obs_data = observable.calculate()

# fetch fluid and slip velocities
boundary_mask = np.squeeze(lb_fluid[:, :, 0].is_boundary.astype(bool))
quivers_boundary = []
quivers_fluid = []
for i, j in itertools.product(range(boundary_mask.shape[0]),
                              range(boundary_mask.shape[1])):
    v_fluid = lb_fluid[i, j, 0].velocity
    if boundary_mask[i, j]:
        quivers_boundary.append([i, j, v_fluid[0], v_fluid[1]])
    else:
        quivers_fluid.append([i, j, v_fluid[0], v_fluid[1]])

# prepare canvas
plt.rcParams.update({'font.size': 16})
fig1 = plt.figure()
fig2 = plt.figure()
ax1 = fig1.add_subplot(111)
ax2 = fig2.add_subplot(111)

# plot velocity as a function of distance
profile_r = observable.bin_centers().reshape([-1, 3])[:, 0]
profile_v = (obs_data - obs_data_baseline).reshape([-1, 3])
ax1.plot(profile_r, profile_v[:, 1])
y_formatter = matplotlib.ticker.ScalarFormatter()
y_formatter.set_powerlimits((-1e-2, 1e-2))
ax1.yaxis.set_major_formatter(y_formatter)
ax1.set(xlabel='Distance from cylinder center', ylabel='Fluid velocity')

# plot boundary geometry
cmap = matplotlib.colors.ListedColormap(['white', 'silver', 'silver'])
cmap_bounds = [0, 1, 2]
cmap_norm = matplotlib.colors.BoundaryNorm(cmap_bounds, cmap.N)
ax2.imshow(boundary_mask.T, origin='lower', interpolation='nearest', cmap=cmap,
           norm=cmap_norm)

# add grid lines based on minor ticks
minor_locator = matplotlib.ticker.FixedLocator(np.arange(0.5, grid_size[0], 1))
ax2.xaxis.set_minor_locator(minor_locator)
ax2.yaxis.set_minor_locator(minor_locator)
ax2.tick_params(axis='both', which='minor', length=0)
ax2.grid(which='minor', color='w', linestyle='-', linewidth=1.2, zorder=2)

# remove major ticks
ax2.set_xticks([])
ax2.set_yticks([])

# add cylinder radii
# circle_in = plt.Circle(
#    cyl_center[:2] / agrid - agrid, cylinder_in.radius / agrid,
#    color='r', fill=False, zorder=3)
# circle_out = plt.Circle(
#    cyl_center[:2] / agrid - agrid, cylinder_out.radius / agrid,
#    color='r', fill=False, zorder=3)
# ax2.add_patch(circle_in)
# ax2.add_patch(circle_out)

# plot velocity field
quivers_boundary = np.array(quivers_boundary)
quivers_fluid = np.array(quivers_fluid)
ax2.quiver(quivers_boundary[:, 0], quivers_boundary[:, 1], quivers_boundary[:, 2],
           quivers_boundary[:, 3], scale=.25, width=0.003, color='black',
           zorder=4, label='slip velocity')
ax2.quiver(quivers_fluid[:, 0], quivers_fluid[:, 1], quivers_fluid[:, 2],
           quivers_fluid[:, 3], scale=.25, width=0.003, color='royalblue',
           zorder=4, label='fluid velocity')
ax2.set(xlabel='x-axis', ylabel='y-axis')
ax2.legend(framealpha=1, loc='upper right')

plt.tight_layout()
plt.show()
