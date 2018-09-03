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
import unittest as ut
import sys
import math
import numpy as np

import espressomd
import espressomd.electrokinetics
import espressomd.shapes
import ek_common
from tests_common import DynamicDict

##########################################################################
#                              Set up the System                               #
##########################################################################
# Set the slit pore geometry the width is the non-periodic part of the geometry
# the padding is used to ensure that there is no field inside outside the slit

params_x = DynamicDict([
            ('box_x', 6.0),
            ('box_y', 6.0),
            ('width', 50.0),
            ('padding', 6.0),
            ('box_z', 'width + 2 * padding'),
            ('agrid', 1.0),
            ('dt', 1.0/7),
            ('force', 0.13),
            ('sigma', -0.05),
            ('viscosity_kinematic', 2.3),
            ('friction', 4.3),
            ('temperature', 2.9),
            ('bjerrum_length', 0.47),
            ('integration_length', 30000),
            ('density_water', 26.15),
            ('density_counterions', '-2.0 * sigma / width'),
            ('valency', 1.0),
            ('ext_force_density', '[force, 0.0, 0.0]'),
            ('wall_normal_1', [0, 0, 1]),
            ('wall_normal_2', [0, 0, -1]),
            ('periodic_length',  'box_z - padding')
])

def bisection(params) :
    # initial parameters for bisection scheme
    size = math.pi / (2.0 * params['width'])
    pnt0 = 0.0
    pntm = pnt0 + size
    pnt1 = pnt0 + 1.9 * size

    # the bisection scheme
    tol = 1.0e-08
    while (size > tol):
        val0 = ek_common.solve(pnt0, params['width'], params['bjerrum_length'], params['sigma'], params['valency'])
        val1 = ek_common.solve(pnt1, params['width'], params['bjerrum_length'], params['sigma'], params['valency'])
        valm = ek_common.solve(pntm, params['width'], params['bjerrum_length'], params['sigma'], params['valency'])

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
    return pntm

@ut.skipIf(not espressomd.has_features(["ELECTROKINETICS", "EK_BOUNDARIES"]),
           "Features not available, skipping test!")
class ek_eof_one_species_x(ut.TestCase):
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    system.seed = system.cell_system.get_state()['n_nodes'] * [1234]

    def tearDown(self):
        self.system.actors.clear()

    def test_x(self):
        self.run_test_x(params_x)

    def run_test_x(self, params):
        system = self.system
        system.box_l = [params['box_x'], params['box_y'], params['box_z']]
        system.time_step = params['dt']
        system.cell_system.skin = 0.1
        system.thermostat.turn_off()

        # Set up the (LB) electrokinetics fluid
        ek = espressomd.electrokinetics.Electrokinetics(
            agrid=params['agrid'],
            lb_density=params['density_water'],
            viscosity=params['viscosity_kinematic'],
            friction=params['friction'],
            T=params['temperature'],
            prefactor=params['bjerrum_length'] * params['temperature'],
            stencil="linkcentered")

        counterions = espressomd.electrokinetics.Species(
            density=params['density_counterions'],
            D=0.3,
            valency=params['valency'],
            ext_force_density=params['ext_force_density'])
        ek.add_species(counterions)

        # Set up the walls confining the fluid and carrying charge
        ek_wall1 = espressomd.ekboundaries.EKBoundary(
            charge_density=params['sigma'] / (params['agrid'] * params['padding']), shape=espressomd.shapes.Wall(normal=params['wall_normal_1'], dist=params['padding']))
        ek_wall2 = espressomd.ekboundaries.EKBoundary(charge_density=params['sigma'] / (
            params['agrid'] * params['padding']), shape=espressomd.shapes.Wall(normal=params['wall_normal_2'], dist=-(params['padding'] + params['width'])))
        system.ekboundaries.add(ek_wall1)
        system.ekboundaries.add(ek_wall2)
        system.actors.add(ek)

        # Integrate the system
        system.integrator.run(params['integration_length'])

        # compare the various quantities to the analytic results
        total_velocity_difference = 0.0
        total_density_difference = 0.0
        total_pressure_difference_xx = 0.0
        total_pressure_difference_yy = 0.0
        total_pressure_difference_zz = 0.0
        total_pressure_difference_xy = 0.0
        total_pressure_difference_yz = 0.0
        total_pressure_difference_xz = 0.0

        # obtain the desired xi value
        xi = bisection(params)
        for i in range(int(params['box_z'] / params['agrid'])):
            if (i * params['agrid'] >= params['padding'] and i * params['agrid'] < params['periodic_length']):
                xvalue = i * params['agrid'] - params['padding']
                position = i * params['agrid'] - params['padding'] - params['width'] / 2.0 + params['agrid'] / 2.0

                # density
                measured_density = counterions[int(
                    params['box_x'] / (2 * params['agrid'])), int(params['box_y'] / (2 * params['agrid'])), i].density
                calculated_density = ek_common.density(position, xi, params['bjerrum_length'])
                density_difference = abs(measured_density - calculated_density)
                total_density_difference = total_density_difference + \
                    density_difference

                # velocity
                measured_velocity = ek[int(
                    params['box_x'] / (2 * params['agrid'])), int(params['box_y'] / (2 * params['agrid'])), i].velocity[0]
                calculated_velocity = ek_common.velocity(
                    position,
                    xi,
                    params['width'],
                    params['bjerrum_length'],
                    params['force'],
                    params['viscosity_kinematic'],
                    params['density_water'])
                velocity_difference = abs(
                    measured_velocity - calculated_velocity)
                total_velocity_difference = total_velocity_difference + \
                    velocity_difference

                # diagonal pressure tensor
                measured_pressure_xx = ek[int(
                    params['box_x'] / (2 * params['agrid'])), int(params['box_y'] / (2 * params['agrid'])), i].pressure[(0, 0)]
                calculated_pressure_xx = ek_common.hydrostatic_pressure(
                    ek, position, xi, params['bjerrum_length'], (0, 0), params['box_x'], params['box_y'], params['box_z'], params['agrid'])
                measured_pressure_yy = ek[int(
                    params['box_x'] / (2 * params['agrid'])), int(params['box_y'] / (2 * params['agrid'])), i].pressure[(1, 1)]
                calculated_pressure_yy = ek_common.hydrostatic_pressure(
                    ek, position, xi, params['bjerrum_length'], (1, 1), params['box_x'], params['box_y'], params['box_z'], params['agrid'])
                measured_pressure_zz = ek[int(
                    params['box_x'] / (2 * params['agrid'])), int(params['box_y'] / (2 * params['agrid'])), i].pressure[(2, 2)]
                calculated_pressure_zz = ek_common.hydrostatic_pressure(
                    ek, position, xi, params['bjerrum_length'], (2, 2), params['box_x'], params['box_y'], params['box_z'], params['agrid'])

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
                    params['box_x'] / (2 * params['agrid'])), int(params['box_y'] / (2 * params['agrid'])), i].pressure[(0, 1)]
                calculated_pressure_xy = 0.0
                pressure_difference_xy = abs(
                    measured_pressure_xy - calculated_pressure_xy)
                total_pressure_difference_xy = total_pressure_difference_xy + \
                    pressure_difference_xy

                # yz component pressure tensor
                measured_pressure_yz = ek[int(
                    params['box_x'] / (2 * params['agrid'])), int(params['box_y'] / (2 * params['agrid'])), i].pressure[(1, 2)]
                calculated_pressure_yz = 0.0
                pressure_difference_yz = abs(
                    measured_pressure_yz - calculated_pressure_yz)
                total_pressure_difference_yz = total_pressure_difference_yz + \
                    pressure_difference_yz

            # xz component pressure tensor
                measured_pressure_xz = ek[int(
                    params['box_x'] / (2 * params['agrid'])), int(params['box_y'] / (2 * params['agrid'])), i].pressure[(0, 2)]
                calculated_pressure_xz = ek_common.pressure_tensor_offdiagonal(
                    position, xi, params['bjerrum_length'], params['force'])
                pressure_difference_xz = abs(
                    measured_pressure_xz - calculated_pressure_xz)
                total_pressure_difference_xz = total_pressure_difference_xz + \
                    pressure_difference_xz

        scale_factor = params['agrid'] / params['width']
        total_density_difference *= scale_factor
        total_velocity_difference *= scale_factor
        total_pressure_difference_xx *= scale_factor
        total_pressure_difference_yy *= scale_factor
        total_pressure_difference_zz *= scale_factor
        total_pressure_difference_xy *= scale_factor
        total_pressure_difference_yz *= scale_factor
        total_pressure_difference_xz *= scale_factor

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
