#
# Copyright (C) 2010-2023 The ESPResSo project
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

import unittest as ut
import unittest_decorators as utx
import espressomd.electrostatics

import numpy as np
import itertools

TIME_STEP = 1e-100


class ElcTest:
    system = espressomd.System(box_l=[1.] * 3, time_step=TIME_STEP)
    system.cell_system.skin = 0.0

    def tearDown(self):
        self.system.part.clear()
        self.system.electrostatics.clear()

    def test_finite_potential_drop(self):
        system = self.system

        GAP = np.array([0., 0., 3.])
        BOX_L = np.array(3 * [10.]) + GAP
        POTENTIAL_DIFFERENCE = -3.

        system.box_l = BOX_L

        p1 = system.part.add(pos=[0, 0, 1], q=+1)
        p2 = system.part.add(pos=[0, 0, 9], q=-1)

        p3m = self.p3m_class(
            # zero is not allowed
            prefactor=1e-100,
            mesh=32,
            cao=5,
            accuracy=1e-3,
        )
        elc = espressomd.electrostatics.ELC(
            actor=p3m,
            gap_size=GAP[2],
            maxPWerror=1e-3,
            delta_mid_top=-1,
            delta_mid_bot=-1,
            const_pot=True,
            pot_diff=POTENTIAL_DIFFERENCE,
        )

        system.electrostatics.solver = elc

        # Calculated energy
        U_elc = system.analysis.energy()['coulomb']

        # Expected E-Field is voltage drop over the box
        E_expected = POTENTIAL_DIFFERENCE / (BOX_L[2] - GAP[2])
        # Expected potential is -E_expected * z, so
        U_expected = -E_expected * (p1.pos[2] * p1.q + p2.pos[2] * p2.q)

        system.integrator.run(0)

        self.assertAlmostEqual(U_elc, U_expected)
        self.assertAlmostEqual(p1.f[2] / p1.q, E_expected)
        self.assertAlmostEqual(p2.f[2] / p2.q, E_expected)

        # Check if error is thrown when particles enter the ELC gap
        # positive direction
        p1.pos = [BOX_L[0] / 2, BOX_L[1] / 2, BOX_L[2] - GAP[2] / 2]
        with self.assertRaises(Exception):
            self.system.analysis.energy()
        with self.assertRaisesRegex(Exception, 'entered ELC gap region'):
            self.system.integrator.run(2)
        # negative direction
        p1.pos = [BOX_L[0] / 2, BOX_L[1] / 2, -GAP[2] / 2]
        with self.assertRaises(Exception):
            self.system.analysis.energy()
        with self.assertRaisesRegex(Exception, 'entered ELC gap region'):
            self.system.integrator.run(2)

    def test_elc_p3m_madelung(self):
        system = self.system

        n_pairs = 6

        BOX_L = 2 * n_pairs
        ELC_GAP = 0.4 * BOX_L
        system.box_l = [BOX_L, BOX_L, BOX_L + ELC_GAP]

        for j, k, l in itertools.product(range(2 * n_pairs), repeat=3):
            system.part.add(pos=[j + 0.5, k + 0.5, l + 0.5],
                            q=(-1)**(j + k + l))

        p3m = self.p3m_class(
            prefactor=1,
            mesh=[60, 60, 84],
            cao=7,
            r_cut=3.314,
            alpha=1.18,
            accuracy=3e-7,
            tune=False,
        )
        elc = espressomd.electrostatics.ELC(
            actor=p3m,
            gap_size=ELC_GAP,
            maxPWerror=3e-7,
            delta_mid_top=-1,
            delta_mid_bot=-1,
            const_pot=True,
            pot_diff=0,
        )

        system.electrostatics.solver = elc

        MADELUNG = -1.74756459463318219
        U_expected = MADELUNG

        U_elc = 2. * \
            system.analysis.energy()['coulomb'] / len(system.part.all())

        np.testing.assert_allclose(U_elc, U_expected, atol=0., rtol=1e-6)


@utx.skipIfMissingFeatures(["P3M"])
class ElcTestCPU(ElcTest, ut.TestCase):

    p3m_class = espressomd.electrostatics.P3M
    rtol = 1e-7


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["P3M"])
class ElcTestGPU(ElcTest, ut.TestCase):

    p3m_class = espressomd.electrostatics.P3MGPU
    rtol = 4e-6


if __name__ == "__main__":
    ut.main()
