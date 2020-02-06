# Copyright(C) 2011-2019 The ESPResSo project
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

import unittest as ut
import unittest_decorators as utx
import espressomd
from espressomd import electrokinetics
import math

##########################################################################
# Set up the System #
##########################################################################
# Build plates using two ek species.


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["ELECTROKINETICS"])
class ek_charged_plate(ut.TestCase):

    es = espressomd.System(box_l=[1.0, 1.0, 1.0])

    def test(self):
        system = self.es

        # Set parameters
        box_x = 20
        box_y = 20
        box_z = 20
        system.box_l = [box_x, box_y, box_z]
        system.cell_system.skin = 0.2
        system.time_step = 0.1
        system.periodicity = [1, 1, 1]
        bjerrum_length = 2.13569
        agrid = 0.5

        system.thermostat.turn_off()

        # Setup the Fluid
        ek = electrokinetics.Electrokinetics(
            agrid=agrid,
            lb_density=1.0,
            viscosity=1.0,
            friction=1.0,
            T=1.0,
            prefactor=bjerrum_length,
            stencil="linkcentered",
            advection=False,
            es_coupling=True)

        positive_ions = electrokinetics.Species(
            density=0.0, D=0.0, valency=1.0)
        negative_ions = electrokinetics.Species(
            density=0.0, D=0.0, valency=-1.0)
        ek.add_species(positive_ions)
        ek.add_species(negative_ions)
        system.actors.add(ek)

        ##################################################################
        # X
        # Setup EK species
        for i in range(int(box_y / agrid)):
            for j in range(int(box_z / agrid)):
                positive_ions[10, i, j].density = 1.0 / agrid
                negative_ions[30, i, j].density = 1.0 / agrid

        # Setup MD particle and integrate
        system.part.add(id=0, pos=[0, 0, 0], q=-1.0, type=0)
        force_difference = 0.0

        for i in range(7, 14):
            system.part[0].pos = [i, 0, 0]
            system.integrator.run(0)

            # Check Force
            expected_force = -2 * math.pi * bjerrum_length
            particle_force = system.part[0].f
            if abs(expected_force - particle_force[0]) > force_difference:
                force_difference = abs(expected_force - particle_force[0])

        self.assertLess(force_difference, 1.0e-04,
                        "Force accuracy in X not achieved, allowed deviation: "
                        "1.0e-04, measured: {}".format(force_difference))

        # Unset species
        for i in range(int(box_y / agrid)):
            for j in range(int(box_z / agrid)):
                positive_ions[10, i, j].density = 0.0
                negative_ions[30, i, j].density = 0.0

        ##################################################################
        # Y
        # Setup EK species
        for i in range(int(box_x / agrid)):
            for j in range(int(box_z / agrid)):
                positive_ions[i, 10, j].density = 1.0 / agrid
                negative_ions[i, 30, j].density = 1.0 / agrid

        # Setup MD particle and integrate
        force_difference = 0.0

        for i in range(7, 14):
            system.part[0].pos = [0, i, 0]
            system.integrator.run(0)

            # Check Force
            expected_force = -2 * math.pi * bjerrum_length
            particle_force = system.part[0].f
            if abs(expected_force - particle_force[1]) > force_difference:
                force_difference = abs(expected_force - particle_force[1])

        self.assertLess(force_difference, 1.0e-04,
                        "Force accuracy in Y not achieved, allowed deviation: "
                        "1.0e-04, measured: {}".format(force_difference))

        # Unset species
        for i in range(int(box_x / agrid)):
            for j in range(int(box_z / agrid)):
                positive_ions[i, 10, j].density = 0.0
                negative_ions[i, 30, j].density = 0.0

        ##################################################################
        # Y
        # Setup EK species
        for i in range(int(box_x / agrid)):
            for j in range(int(box_y / agrid)):
                positive_ions[i, j, 10].density = 1.0 / agrid
                negative_ions[i, j, 30].density = 1.0 / agrid

        # Setup MD particle and integrate
        force_difference = 0.0

        for i in range(7, 14):
            system.part[0].pos = [0, 0, i]
            system.integrator.run(0)

            # Check Force
            expected_force = -2 * math.pi * bjerrum_length
            particle_force = system.part[0].f
            if abs(expected_force - particle_force[2]) > force_difference:
                force_difference = abs(expected_force - particle_force[2])

        self.assertLess(force_difference, 1.0e-04,
                        "Force accuracy in Z not achieved, allowed deviation: "
                        "1.0e-04, measured: {}".format(force_difference))

        # Unset species
        for i in range(int(box_x / agrid)):
            for j in range(int(box_y / agrid)):
                positive_ions[i, j, 10].density = 0.0
                negative_ions[i, j, 30].density = 0.0


if __name__ == "__main__":
    ut.main()
