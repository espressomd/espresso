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
import unittest as ut
import unittest_decorators as utx
import espressomd


@utx.skipIfMissingFeatures(["P3M", "EXTERNAL_FORCES"])
class test_icc(ut.TestCase):

    def runTest(self):
        from espressomd.electrostatics import P3M
        from espressomd.electrostatic_extensions import ICC

        S = espressomd.System(box_l=[1.0, 1.0, 1.0])
        # Parameters
        box_l = 20.0
        nicc = 10
        q_test = 10.0
        q_dist = 5.0

        # System
        S.box_l = [box_l, box_l, box_l + 5.0]
        S.cell_system.skin = 0.4
        S.time_step = 0.01

        # ICC particles
        nicc_per_electrode = nicc * nicc
        nicc_tot = 2 * nicc_per_electrode
        iccArea = box_l * box_l / nicc_per_electrode

        iccNormals = []
        iccAreas = []
        iccSigmas = []
        iccEpsilons = []

        l = box_l / nicc
        for xi in range(nicc):
            for yi in range(nicc):
                S.part.add(pos=[l * xi, l * yi, 0], q=-0.0001, fix=[1, 1, 1])
                iccNormals.append([0, 0, 1])

        for xi in range(nicc):
            for yi in range(nicc):
                S.part.add(pos=[l * xi, l * yi, box_l],
                           q=0.0001, fix=[1, 1, 1])
                iccNormals.append([0, 0, -1])

        iccAreas.extend([iccArea] * nicc_tot)
        iccSigmas.extend([0] * nicc_tot)
        iccEpsilons.extend([10000000] * nicc_tot)

        # Test Dipole
        b2 = box_l * 0.5
        S.part.add(pos=[b2, b2, b2 - q_dist / 2], q=q_test, fix=[1, 1, 1])
        S.part.add(pos=[b2, b2, b2 + q_dist / 2], q=-q_test, fix=[1, 1, 1])

        # Actors
        p3m = P3M(prefactor=1, mesh=32, cao=7, accuracy=1e-5)
        icc = ICC(
            n_icc=nicc_tot,
            convergence=1e-6,
            relaxation=0.75,
            ext_field=[0, 0, 0],
            max_iterations=100,
            first_id=0,
            eps_out=1,
            normals=iccNormals,
            areas=iccAreas,
            sigmas=iccSigmas,
            epsilons=iccEpsilons)

        S.actors.add(p3m)
        S.actors.add(icc)

        # Run
        S.integrator.run(0)

        # Analyze
        QL = sum(S.part[:nicc_per_electrode].q)
        QR = sum(S.part[nicc_per_electrode:nicc_tot].q)

        testcharge_dipole = q_test * q_dist
        induced_dipole = 0.5 * (abs(QL) + abs(QR)) * box_l

        # Result
        self.assertAlmostEqual(1, induced_dipole / testcharge_dipole, places=4)

        # Test applying changes
        enegry_pre_change = S.analysis.energy()['total']
        pressure_pre_change = S.analysis.pressure()['total']
        icc.set_params(sigmas=[2.0] * nicc_tot)
        icc.set_params(epsilons=[20.0] * nicc_tot)
        enegry_post_change = S.analysis.energy()['total']
        pressure_post_change = S.analysis.pressure()['total']
        self.assertNotAlmostEqual(enegry_pre_change, enegry_post_change)
        self.assertNotAlmostEqual(pressure_pre_change, pressure_post_change)


if __name__ == "__main__":
    ut.main()
