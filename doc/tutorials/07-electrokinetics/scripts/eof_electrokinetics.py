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
# Initializing espresso modules and the numpy package
import sys
import numpy as np
import espressomd
espressomd.assert_features(["ELECTROKINETICS"])
from espressomd import electrokinetics, shapes

# Set the slit pore geometry the width is the non-periodic part of the geometry
# the padding is used to ensure that there is no field outside the slit

box_x = 6
box_y = 6
width = 50

padding = 1
box_z = width + 2 * padding

system = espressomd.System(box_l=[box_x, box_y, box_z])

# Set the electrokinetic parameters

agrid = 1.0
dt = 0.2
kT = 1.0
bjerrum_length = 0.7095
D = 0.006075
valency = 1.0
viscosity_dynamic = 79.53
density_water = 26.15
sigma = -0.05
ext_force_density = 0.1

# Set the simulation parameters

system.time_step = dt
system.cell_system.skin = 0.2
system.thermostat.turn_off()
integration_length = 2000

# Set up the (LB) electrokinetics fluid
viscosity_kinematic = viscosity_dynamic / density_water
ek = electrokinetics.Electrokinetics(agrid=agrid, lb_density=density_water,
                                     viscosity=viscosity_kinematic, friction=1.0,
                                     T=kT, prefactor=bjerrum_length)

# Set up the charged and neutral species
density_counterions = -2.0 * sigma / width
counterions = electrokinetics.Species(density=density_counterions,
                                      D=D, valency=valency,
                                      ext_force_density=[ext_force_density, 0, 0])

ek.add_species(counterions)

# Set up the walls confining the fluid
ek_wall_left = espressomd.ekboundaries.EKBoundary(charge_density=sigma / agrid,
                                                  shape=shapes.Wall(
                                                  normal=[0, 0, 1],
                                                  dist=padding))
ek_wall_right = espressomd.ekboundaries.EKBoundary(
    charge_density=sigma / agrid,
                                           shape=shapes.Wall(normal=[0, 0, -1],
                                                             dist=-(padding + width)))

system.ekboundaries.add(ek_wall_left)
system.ekboundaries.add(ek_wall_right)

system.actors.add(ek)

# Integrate the system
for i in range(100):
    system.integrator.run(integration_length)
    sys.stdout.write("\rintegration: %i%%" % (i + 1))
    sys.stdout.flush()

# Output
position_list = []
density_list = []
velocity_list = []
pressure_xz_list = []

for i in range(int(box_z / agrid)):
    if (i * agrid >= padding) and (i * agrid < box_z - padding):
        position = i * agrid - padding - width / 2.0 + agrid / 2.0
        position_list.append(position)

        # density
        density_list.append(
            counterions[box_x / (2 * agrid), box_y / (2 * agrid), i].density)

        # velocity
        velocity_list.append(
            ek[box_x / (2 * agrid), box_y / (2 * agrid), i].velocity[0])

        # xz component pressure tensor
        pressure_xz_list.append(
            ek[box_x / (2 * agrid), box_y / (2 * agrid), i].pressure[0, 2])

np.savetxt("eof_electrokinetics.dat", np.column_stack((position_list,
                                                       density_list,
                                                       velocity_list,
                                                       pressure_xz_list)),
           header="#position calculated_density calculated_velocity calculated_pressure_xz")
