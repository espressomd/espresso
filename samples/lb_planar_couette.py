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
Simulate the flow profile of a lattice-Boltzmann fluid between two
shear planes with Lees-Edwards boundary conditions and compare it
to the analytical solution.
"""

import espressomd
import espressomd.lb
import espressomd.lees_edwards

import numpy as np
import matplotlib.pyplot as plt

required_features = ["WALBERLA"]
espressomd.assert_features(required_features)


def analytical(x, t, nu, v, h, k_max):
    """
    Analytical solution with Fourier series of the Navier-Stokes equation.

    Parameters
    ----------
    x : :obj:`float`
        Height within the channel
    t : :obj:`float`
        Time since the start up of the shear flow
    nu: :obj:`float`
        Kinematic viscosity
    v: :obj:`float`
        Shearing velocity
    h : :obj:`float`
        Distance between shear planes
    k_max : :obj:`int`
        Upper limit of sums for sinus series
    """
    u = x / h - 0.5
    for k in np.arange(1, k_max + 1):
        wave = 2 * np.pi * k / h
        u += np.exp(-nu * wave ** 2 * t) * np.sin(wave * x) / (np.pi * k)
    return v * u


# LB and LE parameters
nu = 1. / 6.
h = 64.0
v = 0.02
k_max = 100

system = espressomd.System(box_l=[h, 64, 1])
system.time_step = 1.
system.cell_system.skin = 0.1
system.cell_system.set_n_square()

system.lees_edwards.set_boundary_conditions(
    shear_direction="x", shear_plane_normal="y",
    protocol=espressomd.lees_edwards.LinearShear(
        shear_velocity=v, initial_pos_offset=0.0, time_0=0.0))

lbf = espressomd.lb.LBFluidWalberla(
    agrid=1., density=1., kinematic_viscosity=nu, tau=1.)
system.actors.add(lbf)

# sampling
time_breakpoints = [50, 200, 500, 2000]
pos_breakpoints = 256
for steps in time_breakpoints:
    steps -= int(system.time) - 1
    system.integrator.run(steps)
    time = system.time - 1.
    position_ref = np.linspace(0.5, 63.5, pos_breakpoints)
    position_lbf = np.linspace(0.5, 63.5, 64)
    velocity_ref = analytical(position_ref, time, nu, v, h, k_max)
    velocity_lbf = np.copy(lbf[5, :, 0].velocity[:, 0].reshape([-1]))
    ax = plt.gca()
    color = next(ax._get_lines.prop_cycler)['color']
    plt.plot(velocity_ref, position_ref, '-', color=color,
             label=f"Analytical solution at t={time:.0f}")
    plt.plot(velocity_lbf, position_lbf, 'o', color=color,
             label=f"Simulated profile at t={time:.0f}")

plt.xlabel('shear velocity')
plt.ylabel('y-position')
# format legend in 2 columns
ax = plt.gca()
handles, labels = ax.get_legend_handles_labels()
handles, labels = zip(*sorted(zip(handles, labels), key=lambda x: x[1][0]))
ax.legend(handles, labels, ncol=2)
plt.show()
