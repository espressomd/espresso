# Copyright (C) 2011-2018 The ESPResSo project
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
import sys
import math
import unittest as ut
import numpy as np

import espressomd
import espressomd.electrokinetics
import espressomd.shapes
import ek_common

##########################################################################
#                              Set up the System                               #
##########################################################################
# Set the slit pore geometry the width is the non-periodic part of the geometry
# the padding is used to ensure that there is no field inside outside the slit


@ut.skipIf(not espressomd.has_features(["ELECTROKINETICS", "EK_BOUNDARIES"]),
           "Features not available, skipping test!")
class ek_eof_one_species_x(ut.TestCase):
    es = espressomd.System(box_l=[1.0, 1.0, 1.0])
    es.seed = es.cell_system.get_state()['n_nodes'] * [1234]

    def test(self):
        system = self.es

        pi = math.pi
        box_x = 6
        box_y = 6
        width = 20

        padding = 6
        box_z = width + 2 * padding

        # Set the electrokinetic parameters
        agrid = 1.0
        dt = 1.0 / 7.
        force = 0.13
        sigma = -0.05
        viscosity_kinematic = 2.3
        friction = 4.3
        temperature = 2.9
        bjerrum_length = 0.47

        temperature_LB = agrid * agrid / (3.0 * dt * dt)
        kB_LB = 1.0
        cs_squared = (1.0 / 3.0) * (agrid * agrid / (dt * dt))

        system.box_l = [box_x, box_y, box_z]

        # Set the simulation parameters
        system.time_step = dt
        system.cell_system.skin = 0.1
        system.thermostat.turn_off()
        integration_length = 4000

        # Set up the charged and neutral species
        density_water = 26.15
        density_counterions = -2.0 * float(sigma) / float(width)
        valency = 1.0

        # Set up the (LB) electrokinetics fluid
        ek = espressomd.electrokinetics.Electrokinetics(
            agrid=agrid,
            lb_density=density_water,
            viscosity=viscosity_kinematic,
            friction=friction,
            T=temperature,
            prefactor=bjerrum_length * temperature,
            stencil="nonlinear")

        counterions = espressomd.electrokinetics.Species(
            density=density_counterions,
            D=0.3,
            valency=valency,
            ext_force_density=[
                force,
                0,
                0])
        ek.add_species(counterions)

        # Set up the walls confining the fluid and carrying charge
        ek_wall1 = espressomd.ekboundaries.EKBoundary(
            charge_density=sigma / (
                agrid * padding),
            shape=espressomd.shapes.Wall(
                normal=[
                    0,
                    0,
                    1],
                dist=padding))
        system.ekboundaries.add(ek_wall1)
        ek_wall2 = espressomd.ekboundaries.EKBoundary(charge_density=sigma / (
            agrid * padding), shape=espressomd.shapes.Wall(normal=[0, 0, -1], dist=-(padding + width)))
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
            val0 = ek_common.solve(pnt0, width, bjerrum_length, sigma, valency)
            val1 = ek_common.solve(pnt1, width, bjerrum_length, sigma, valency)
            valm = ek_common.solve(pntm, width, bjerrum_length, sigma, valency)

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

        for i in range(int(box_z / agrid)):

            if (i * agrid >= padding and i * agrid < box_z - padding):
                xvalue = i * agrid - padding
                position = i * agrid - padding - width / 2.0 + agrid / 2.0

                # density
                measured_density = counterions[int(
                    box_x / (2 * agrid)), int(box_y / (2 * agrid)), i].density
                calculated_density = ek_common.density(
                    position, xi, bjerrum_length)
                density_difference = abs(measured_density - calculated_density)
                total_density_difference = total_density_difference + \
                    density_difference

                # velocity
                measured_velocity = ek[int(
                    box_x / (2 * agrid)), int(box_y / (2 * agrid)), i].velocity[0]
                calculated_velocity = ek_common.velocity(
                    position,
                    xi,
                    width,
                    bjerrum_length,
                    force,
                    viscosity_kinematic,
                    density_water)
                velocity_difference = abs(
                    measured_velocity - calculated_velocity)
                total_velocity_difference = total_velocity_difference + \
                    velocity_difference

                # diagonal pressure tensor
                measured_pressure_xx = ek[int(
                    box_x / (2 * agrid)), int(box_y / (2 * agrid)), i].pressure[(0, 0)]
                calculated_pressure_xx = ek_common.hydrostatic_pressure_non_lin(
                    ek, position, xi, bjerrum_length, (0, 0), box_x, box_y, box_z, agrid, temperature)
                measured_pressure_yy = ek[int(
                    box_x / (2 * agrid)), int(box_y / (2 * agrid)), i].pressure[(1, 1)]
                calculated_pressure_yy = ek_common.hydrostatic_pressure_non_lin(
                    ek, position, xi, bjerrum_length, (1, 1), box_x, box_y, box_z, agrid, temperature)
                measured_pressure_zz = ek[int(
                    box_x / (2 * agrid)), int(box_y / (2 * agrid)), i].pressure[(2, 2)]
                calculated_pressure_zz = ek_common.hydrostatic_pressure_non_lin(
                    ek, position, xi, bjerrum_length, (2, 2), box_x, box_y, box_z, agrid, temperature)

                pressure_difference_xx = abs(
                    measured_pressure_xx - calculated_pressure_xx)
                pressure_difference_yy = abs(
                    measured_pressure_yy - calculated_pressure_yy)
                pressure_difference_zz = abs(
                    measured_pressure_zz - calculated_pressure_zz)

                total_pressure_difference_xx = total_pressure_difference_xx + \
                    pressure_difference_xx
                total_pressure_difference_yy = total_pressure_difference_yy + \
                    pressure_difference_yy
                total_pressure_difference_zz = total_pressure_difference_zz + \
                    pressure_difference_zz

                # xy component pressure tensor
                measured_pressure_xy = ek[int(
                    box_x / (2 * agrid)), int(box_y / (2 * agrid)), i].pressure[(0, 1)]
                calculated_pressure_xy = 0.0
                pressure_difference_xy = abs(
                    measured_pressure_xy - calculated_pressure_xy)
                total_pressure_difference_xy = total_pressure_difference_xy + \
                    pressure_difference_xy

                # yz component pressure tensor
                measured_pressure_yz = ek[int(
                    box_x / (2 * agrid)), int(box_y / (2 * agrid)), i].pressure[(1, 2)]
                calculated_pressure_yz = 0.0
                pressure_difference_yz = abs(
                    measured_pressure_yz - calculated_pressure_yz)
                total_pressure_difference_yz = total_pressure_difference_yz + \
                    pressure_difference_yz

                # xz component pressure tensor
                measured_pressure_xz = ek[int(
                    box_x / (2 * agrid)), int(box_y / (2 * agrid)), i].pressure[(0, 2)]
                calculated_pressure_xz = ek_common.pressure_tensor_offdiagonal(
                    position, xi, bjerrum_length, force)
                pressure_difference_xz = abs(
                    measured_pressure_xz - calculated_pressure_xz)
                total_pressure_difference_xz = total_pressure_difference_xz + \
                    pressure_difference_xz

        total_density_difference = agrid * total_density_difference / width
        total_velocity_difference = agrid * total_velocity_difference / width
        total_pressure_difference_xx = agrid * \
            total_pressure_difference_xx / width
        total_pressure_difference_yy = agrid * \
            total_pressure_difference_yy / width
        total_pressure_difference_zz = agrid * \
            total_pressure_difference_zz / width
        total_pressure_difference_xy = agrid * \
            total_pressure_difference_xy / width
        total_pressure_difference_yz = agrid * \
            total_pressure_difference_yz / width
        total_pressure_difference_xz = agrid * \
            total_pressure_difference_xz / width

        self.assertLess(total_density_difference, 1.0e-04,
                        "Density accuracy not achieved")
        self.assertLess(total_velocity_difference, 1.0e-04,
                        "Velocity accuracy not achieved")
        self.assertLess(total_pressure_difference_xx, 1.0e-04,
                        "Pressure accuracy xx component not achieved")
        self.assertLess(total_pressure_difference_yy, 1.0e-04,
                        "Pressure accuracy yy component not achieved")
        self.assertLess(total_pressure_difference_zz, 1.0e-04,
                        "Pressure accuracy zz component not achieved")
        self.assertLess(total_pressure_difference_xy, 1.0e-04,
                        "Pressure accuracy xy component not achieved")
        self.assertLess(total_pressure_difference_yz, 1.0e-04,
                        "Pressure accuracy yz component not achieved")
        self.assertLess(total_pressure_difference_xz, 1.0e-04,
                        "Pressure accuracy xz component not achieved")


if __name__ == "__main__":
    ut.main()
