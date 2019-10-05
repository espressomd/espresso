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
import math

# Some common functions used by the ek tests

# root finding function


def solve(xi, d, bjerrum_length, sigma, valency):
    pi = math.pi
    el_char = 1.0
    return xi * math.tan(xi * d / 2.0) + 2.0 * pi * \
        bjerrum_length * sigma / (valency * el_char)

# function to calculate the density


def density(x, xi, bjerrum_length):
    pi = math.pi
    kb = 1.0
    return (xi * xi) / (2.0 * pi * bjerrum_length *
                        math.cos(xi * x) * math.cos(xi * x))

# function to calculate the velocity


def velocity(
        x,
        xi,
        d,
        bjerrum_length,
        force,
        viscosity_kinematic,
        density_water):
    pi = math.pi
    return force * math.log(math.cos(xi * x) / math.cos(xi * d / 2.0)) / \
        (2.0 * pi * bjerrum_length * viscosity_kinematic * density_water)

# function to calculate the nonzero component of the pressure tensor


def pressure_tensor_offdiagonal(x, xi, bjerrum_length, force):
    pi = math.pi
    return force * xi * math.tan(xi * x) / (2.0 * pi * bjerrum_length)

# function to calculate the hydrostatic pressure

# Technically, the LB simulates a compressible fluid, whiches pressure
# tensor contains an additional term on the diagonal, proportional to
# the divergence of the velocity. We neglect this contribution, which
# creates a small error in the direction normal to the wall, which
# should decay with the simulation time.


def hydrostatic_pressure(
        ek,
        x,
        xi,
        bjerrum_length,
        tensor_entry,
        box_x,
        box_y,
        box_z,
        agrid):
    offset = ek[int(box_x / (2 * agrid)), int(box_y / (2 * agrid)),
                int(box_z / (2 * agrid))].pressure[tensor_entry]
    return 0.0 + offset


# variant from the nonlinear tests
def hydrostatic_pressure_non_lin(
        ek,
        x,
        xi,
        bjerrum_length,
        tensor_entry,
        box_x,
        box_y,
        box_z,
        agrid,
        temperature):
    offset = ek[int(box_x / (2 * agrid)), int(box_y / (2 * agrid)),
                int(box_z / (2 * agrid))].pressure[tensor_entry]
    return temperature * xi * xi * \
        math.tan(xi * x) * math.tan(xi * x) / \
        (2.0 * math.pi * bjerrum_length) + offset
