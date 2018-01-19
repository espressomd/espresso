# Copyright (C) 2011,2012,2013,2014,2015,2016 The ESPResSo project
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

from __future__ import print_function
import unittest as ut
import espressomd
import espressomd.electrokinetics
import espressomd.shapes
from espressomd import *
import numpy as np
import sys
import math

##########################################################################
#                              Set up the System                               #
##########################################################################
# Set the slit pore geometry the width is the non-periodic part of the geometry
# the padding is used to ensure that there is no field inside outside the slit

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
        agrid,
        temperature):
    offset = ek[int(box_x / (2 * agrid)), int(box_y / (2 * agrid)),
                int(box_z / (2 * agrid))].pressure[tensor_entry]
    return temperature * xi * xi * \
        math.tan(xi * x) * math.tan(xi * x) / (2.0 * math.pi * bjerrum_length) + offset


@ut.skipIf(not espressomd.has_features(["ELECTROKINETICS", "EK_BOUNDARIES"]),
           "Features not available, skipping test!")
class ek_eof_one_species_x(ut.TestCase):

    es = espressomd.System()

    def test(self):
        system = self.es

        pi = math.pi
        box_z = 6
        box_x = 6
        width = 50

        padding = 6
        box_y = width + 2 * padding

# Set the electrokinetic parameters
        agrid = 1.0 / 3.0
        dt = 1.0 / 13.0
        force = 0.07
        sigma = -0.04
        viscosity_kinematic = 1.7
        friction = 1.9
        temperature = 1.1
        bjerrum_length = 0.8

        temperature_LB = agrid * agrid / (3.0 * dt * dt)
        kB_LB = 1.0
        cs_squared = (1.0 / 3.0) * (agrid * agrid / (dt * dt))

        system.box_l = [box_x, box_y, box_z]


# Set the simulation parameters

        system.time_step = dt
        system.cell_system.skin = 0.1
        system.thermostat.turn_off()
        integration_length = 20000

# Output density, velocity, and pressure tensor profiles

        output_profiles = 1

# Set up the charged and neutral species

        density_water = 26.15
        density_counterions = -2.0 * float(sigma) / float(width)
        valency = 1.0

# Set up the (LB) electrokinetics fluid

        ek = electrokinetics.Electrokinetics(
            agrid=agrid,
            lb_density=density_water,
            viscosity=viscosity_kinematic,
            friction=friction,
            T=temperature,
            prefactor=bjerrum_length*temperature,
            stencil="nonlinear")

        counterions = electrokinetics.Species(
            density=density_counterions,
            D=0.3,
            valency=valency,
            ext_force=[
                0,
                0,
                force])
        ek.add_species(counterions)


# Set up the walls confining the fluid and carrying charge

        ek_wall1 = electrokinetics.EKBoundary(
            charge_density=sigma /
            (padding),
            shape=shapes.Wall(
                normal=[
                    0,
                    1,
                    0],
                dist=padding))
        system.ekboundaries.add(ek_wall1)
        ek_wall2 = electrokinetics.EKBoundary(
            charge_density=sigma / (padding), shape=shapes.Wall(normal=[0, -1, 0], dist=-(padding + width)))
        system.ekboundaries.add(ek_wall2)

        system.actors.add(ek)

# Integrate the system

        system.integrator.run(integration_length)

# compare the various quantities to the analytic results

        total_velocity_difference = 0.0
        total_density_difference = 0.0
        total_pressure_difference_xx = 0.0
        total_pressure_difference_yy = 0.0
        total_pressure_difference_zz = 0.0
        total_pressure_difference_xy = 0.0
        total_pressure_difference_yz = 0.0
        total_pressure_difference_xz = 0.0

# initial parameters for bisection scheme
        size = pi / (2.0 * width)

        pnt0 = 0.0
        pntm = pnt0 + size
        pnt1 = pnt0 + 1.9 * size

# the bisection scheme

        tol = 1.0e-08
        while (size > tol):

            val0 = solve(pnt0, width, bjerrum_length, sigma, valency)
            val1 = solve(pnt1, width, bjerrum_length, sigma, valency)
            valm = solve(pntm, width, bjerrum_length, sigma, valency)

            if (val0 < 0.0 and val1 > 0.0):
                if (valm < 0.0):
                    pnt0 = pntm
                    size = size / 2.0
                    pntm = pnt0 + size
                else:
                    pnt1 = pntm
                    size = size / 2.0
                    pntm = pnt1 - size
            elif (val0 > 0.0 and val1 < 0.0):
                if (valm < 0.0):
                    pnt1 = pntm
                    size = size / 2.0
                    pntm = pnt1 - size
                else:
                    pnt0 = pntm
                    size = size / 2.0
                    pntm = pnt0 + size
            else:
                sys.exit(
                    "Bisection method fails:\nTuning of domain boundaries may be required.")

# obtain the desired xi value

        xi = pntm

        if (output_profiles):
            fp = open("ek_eof_profile.dat", "w")

        for i in range(int(box_y / agrid)):

            if (i * agrid >= padding and i * agrid < box_y - padding):
                xvalue = i * agrid - padding
                position = i * agrid - padding - width / 2.0 + agrid / 2.0

            # density
                measured_density = counterions[int(
                    box_x / (2 * agrid)), i, int(box_z / (2 * agrid))].density
                calculated_density = density(position, xi, bjerrum_length)
                density_difference = abs(measured_density - calculated_density)
                total_density_difference = total_density_difference + density_difference
            # velocity
                measured_velocity = ek[int(
                    box_x / (2 * agrid)), i, int(box_z / (2 * agrid))].velocity[2]
                calculated_velocity = velocity(
                    position,
                    xi,
                    width,
                    bjerrum_length,
                    force,
                    viscosity_kinematic,
                    density_water)
                velocity_difference = abs(
                    measured_velocity - calculated_velocity)
                total_velocity_difference = total_velocity_difference + velocity_difference

            # diagonal pressure tensor

                measured_pressure_xx = ek[int(
                    box_x / (2 * agrid)), i, int(box_z / (2 * agrid))].pressure[(0, 0)]
                calculated_pressure_xx = hydrostatic_pressure(
                    ek, position, xi, bjerrum_length, (0, 0), box_x, box_y, box_z, agrid, temperature)
                measured_pressure_yy = ek[int(
                    box_x / (2 * agrid)), i, int(box_z / (2 * agrid))].pressure[(1, 1)]
                calculated_pressure_yy = hydrostatic_pressure(
                    ek, position, xi, bjerrum_length, (1, 1), box_x, box_y, box_z, agrid, temperature)
                measured_pressure_zz = ek[int(
                    box_x / (2 * agrid)), i, int(box_z / (2 * agrid))].pressure[(2, 2)]
                calculated_pressure_zz = hydrostatic_pressure(
                    ek, position, xi, bjerrum_length, (2, 2), box_x, box_y, box_z, agrid, temperature)

                pressure_difference_xx = abs(
                    measured_pressure_xx - calculated_pressure_xx)
                pressure_difference_yy = abs(
                    measured_pressure_yy - calculated_pressure_yy)
                pressure_difference_zz = abs(
                    measured_pressure_zz - calculated_pressure_zz)

                total_pressure_difference_xx = total_pressure_difference_xx + pressure_difference_xx
                total_pressure_difference_yy = total_pressure_difference_yy + pressure_difference_yy
                total_pressure_difference_zz = total_pressure_difference_zz + pressure_difference_zz

            # xy component pressure tensor
                measured_pressure_xy = ek[int(
                    box_x / (2 * agrid)), i, int(box_z / (2 * agrid))].pressure[(0, 1)]
                calculated_pressure_xy = 0.0
                pressure_difference_xy = abs(
                    measured_pressure_xy - calculated_pressure_xy)
                total_pressure_difference_xy = total_pressure_difference_xy + pressure_difference_xy

            # yz component pressure tensor
                measured_pressure_yz = ek[int(
                    box_x / (2 * agrid)), i, int(box_z / (2 * agrid))].pressure[(1, 2)]
                calculated_pressure_yz = pressure_tensor_offdiagonal(
                    position, xi, bjerrum_length, force)
                pressure_difference_yz = abs(
                    measured_pressure_yz - calculated_pressure_yz)
                total_pressure_difference_yz = total_pressure_difference_yz + pressure_difference_yz

            # xz component pressure tensor
                measured_pressure_xz = ek[int(
                    box_x / (2 * agrid)), i, int(box_z / (2 * agrid))].pressure[(0, 2)]
                calculated_pressure_xz = 0.0
                pressure_difference_xz = abs(
                    measured_pressure_xz - calculated_pressure_xz)
                total_pressure_difference_xz = total_pressure_difference_xz + pressure_difference_xz

                if (output_profiles):
                    fp.write(
                        "{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(
                            position,
                            measured_density,
                            calculated_density,
                            measured_velocity,
                            calculated_velocity,
                            measured_pressure_xy,
                            calculated_pressure_xy,
                            measured_pressure_yz,
                            calculated_pressure_yz,
                            measured_pressure_xz,
                            calculated_pressure_xz,
                            measured_pressure_xx,
                            calculated_pressure_xx,
                            measured_pressure_yy,
                            calculated_pressure_yy,
                            measured_pressure_zz,
                            calculated_pressure_zz))

        if (output_profiles):
            fp.close

        total_density_difference = agrid * total_density_difference / width
        total_velocity_difference = agrid * total_velocity_difference / width
        total_pressure_difference_xx = agrid * total_pressure_difference_xx / width
        total_pressure_difference_yy = agrid * total_pressure_difference_yy / width
        total_pressure_difference_zz = agrid * total_pressure_difference_zz / width
        total_pressure_difference_xy = agrid * total_pressure_difference_xy / width
        total_pressure_difference_yz = agrid * total_pressure_difference_yz / width
        total_pressure_difference_xz = agrid * total_pressure_difference_xz / width

        print("Density deviation: {}".format(total_density_difference))
        print("Velocity deviation: {}".format(total_velocity_difference))
        print("Pressure deviation xx component: {}".format(
            total_pressure_difference_xx))
        print("Pressure deviation yy component: {}".format(
            total_pressure_difference_yy))
        print("Pressure deviation zz component: {}".format(
            total_pressure_difference_zz))
        print("Pressure deviation xy component: {}".format(
            total_pressure_difference_xy))
        print("Pressure deviation yz component: {}".format(
            total_pressure_difference_yz))
        print("Pressure deviation xz component: {}".format(
            total_pressure_difference_xz))

        self.assertLess(total_density_difference, 1.5e-06,
                        "Density accuracy not achieved")
        self.assertLess(total_velocity_difference, 4.0e-07,
                        "Velocity accuracy not achieved")
        self.assertLess(total_pressure_difference_xx, 5.0e-06,
                        "Pressure accuracy xx component not achieved")
        self.assertLess(total_pressure_difference_yy, 5.0e-06,
                        "Pressure accuracy yy component not achieved")
        self.assertLess(total_pressure_difference_zz, 4.0e-06,
                        "Pressure accuracy zz component not achieved")
        self.assertLess(total_pressure_difference_xy, 1.0e-10,
                        "Pressure accuracy xy component not achieved")
        self.assertLess(total_pressure_difference_yz, 1.5e-06,
                        "Pressure accuracy yz component not achieved")
        self.assertLess(total_pressure_difference_xz, 1.0e-10,
                        "Pressure accuracy xz component not achieved")


if __name__ == "__main__":
    ut.main()
