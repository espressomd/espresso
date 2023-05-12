#
# Copyright (C) 2010-2022 The ESPResSo project
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
Simulate the flow of a lattice-Boltzmann fluid past a cylinder,
obtain the velocity profile in polar coordinates and compare it
to the analytical solution.
"""
import matplotlib.pyplot as plt

import espressomd
import espressomd.lb
import espressomd.observables
import espressomd.shapes
import espressomd.accumulators
import espressomd.math

required_features = ["WALBERLA"]
espressomd.assert_features(required_features)

system = espressomd.System(box_l=[10.0, 10.0, 5.0])
system.time_step = 0.01
system.cell_system.skin = 0.4
radius = 4.0
n_steps_warmup = 1000
n_steps = 800

lb_fluid = espressomd.lb.LBFluidWalberla(
    agrid=1.0, density=1.0, kinematic_viscosity=1.0, tau=0.01,
    ext_force_density=[0, 0, 0.15], kT=0.0)
system.actors.add(lb_fluid)
system.thermostat.set_lb(LB_fluid=lb_fluid, seed=23)
ctp = espressomd.math.CylindricalTransformationParameters(
    center=[5.0, 5.0, 0.0],
    axis=[0, 0, 1],
    orientation=[1, 0, 0])
fluid_obs = espressomd.observables.CylindricalLBVelocityProfile(
    transform_params=ctp,
    n_r_bins=16,
    max_r=radius,
    min_z=0.0,
    max_z=system.box_l[2] / lb_fluid.agrid,
    sampling_density=0.1)
cylinder_shape = espressomd.shapes.Cylinder(
    center=[5.0, 5.0, 5.0],
    axis=[0, 0, 1],
    direction=-1,
    radius=radius,
    length=20.0)
lb_fluid.add_boundary_from_shape(cylinder_shape)

# equilibrate fluid
system.integrator.run(n_steps_warmup)
# sample
accumulator = espressomd.accumulators.MeanVarianceCalculator(obs=fluid_obs)
system.auto_update_accumulators.add(accumulator)
system.integrator.run(n_steps)

lb_fluid_profile = accumulator.mean()


def poiseuille_flow(r, R, ext_force_density):
    return ext_force_density * 1. / 4 * (R**2.0 - r**2.0)


r = fluid_obs.bin_centers()[:, 0, 0, 0]
expected_profile = poiseuille_flow(r, radius, lb_fluid.ext_force_density[2])
plt.plot(r, expected_profile, label='analytical solution')
plt.plot(r, lb_fluid_profile[:, 0, 0, 2], '+', label='LB profile')
plt.xlabel('Distance from pipe center')
plt.ylabel('Fluid velocity')
plt.legend()
plt.show()
