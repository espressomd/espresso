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
# Set the slit pore geometry the width is the non-periodic part of the geometry
# the padding is used to ensure that there is no field outside the slit
import numpy as np

box_x = 6
box_y = 6
width = 50

padding = 6
box_z = width + 2 * padding

# Set the electrokinetic parameters

agrid = 1.0
temperature = 1.0
bjerrum_length = 0.7095
valency = 1.0
viscosity_dynamic = 79.53
density_water = 26.15
sigma = -0.05
force = 0.1

viscosity_kinematic = viscosity_dynamic / density_water
density_counterions = -2.0 * sigma / width

# Calculate the inverse length xi, which is a combination of various
# constants (xi = zeC/2kBT), with C a constant that needs to be
# solved for, or equivalently, xi needs to be solved for

# root finding function


def solve(xi=None, d=None, bjerrum_length=None, sigma=None, valency=None):
    el_char = 1.0
    return xi * np.tan(xi * d / 2.0) + 2.0 * np.pi * \
        bjerrum_length * sigma / (valency * el_char)


size = np.pi / (2.0 * width)

pnt0 = 0.0
pntm = pnt0 + size
pnt1 = pnt0 + 1.9 * size

# the bisection scheme
tol = 1e-8

while size > tol:
    val0 = solve(xi=pnt0, d=width,
                 bjerrum_length=bjerrum_length, sigma=sigma, valency=valency)
    val1 = solve(xi=pnt1, d=width,
                 bjerrum_length=bjerrum_length, sigma=sigma, valency=valency)
    valm = solve(xi=pntm, d=width,
                 bjerrum_length=bjerrum_length, sigma=sigma, valency=valency)

    if (val0 < 0.0) and (val1 > 0.0):
        if valm < 0.0:
            pnt0 = pntm
            size = size / 2.0
            pntm = pnt0 + size
        else:
            pnt1 = pntm
            size = size / 2.0
            pntm = pnt1 - size
    elif (val0 > 0.0) and (val1 < 0.0):
        if valm < 0.0:
            pnt1 = pntm
            size = size / 2.0
            pntm = pnt0 - size
        else:
            pnt0 = pntm
            size = size / 2.0
            pntm = pnt0 - size
    else:
        raise Exception(
            "Bisection method fails:\nTuning of domain boundaries may be required.")


# obtain the desired xi value
xi = pntm

# function to calculate the density


def density(x=None, xi=None, bjerrum_length=None):
    kb = 1.0
    return (xi**2) / (2.0 * np.pi * bjerrum_length * np.cos(xi * x)**2)

# function to calculate the velocity


def velocity(x=None, xi=None, d=None, bjerrum_length=None, force=None,
             viscosity_kinematic=None, density_water=None):
    return force * np.log(np.cos(xi * x) / np.cos(xi * d / 2.0)) / \
        (2.0 * np.pi * bjerrum_length * viscosity_kinematic * density_water)

# function to calculate the nonzero component of the pressure tensor


def pressure_tensor_offdiagonal(x=None, xi=None, bjerrum_length=None,
                                force=None):
    return force * xi * np.tan(xi * x) / (2.0 * np.pi * bjerrum_length)

# function to calculate the hydrostatic pressure

# Technically, the LB simulates a compressible fluid, whose pressure
# tensor contains an additional term on the diagonal, proportional to
# the divergence of the velocity. We neglect this (small) contribution.
# The LB pressure tensor also contains the convective acceleration, which
# we neglect here.


def hydrostatic_pressure(x=None, xi=None, bjerrum_length=None,
                         tensor_entry=None):
    return 0.0


position_list = []
density_list = []
velocity_list = []
pressure_xy_list = []

for i in range(int(box_z / agrid)):
    if (i * agrid >= padding) and (i * agrid < box_z - padding):
        position = i * agrid - padding - width / 2.0 + agrid / 2.0
        position_list.append(position)

        # density
        density_list.append(density(x=position, xi=xi,
                                    bjerrum_length=bjerrum_length))

        # velocity
        velocity_list.append(velocity(x=position, xi=xi,
                                      d=width, bjerrum_length=bjerrum_length,
                                      force=force, viscosity_kinematic=viscosity_kinematic,
                                      density_water=density_water))
        # xz component pressure tensor
        pressure_xy_list.append(pressure_tensor_offdiagonal(x=position, xi=xi,
                                                            bjerrum_length=bjerrum_length, force=force))

np.savetxt(
    "eof_analytical.dat", np.column_stack(
        (position_list, density_list, velocity_list, pressure_xy_list)),
    header="#position calculated_density calculated_velocity calculated_pressure_xy")
