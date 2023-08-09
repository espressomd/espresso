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

import unittest as ut
import unittest_decorators as utx
import espressomd
import numpy as np
import math
import espressomd.interactions
import espressomd.electrostatics

COULOMB_PREFACTOR = 4.01


@utx.skipIfMissingFeatures(["THOLE", "EXTERNAL_FORCES", "P3M"])
class TestThole(ut.TestCase):

    """
    This testcase takes a large box to minimize periodic effects and tests the
    Thole damping non-bonded interaction forces against the analytical result.
    """

    q1 = 1.01
    q2 = -1.01
    thole_s = 1.99

    box_l = 500.0
    system = espressomd.System(box_l=[box_l] * 3)

    def setUp(self):
        self.system.time_step = 0.01
        self.system.cell_system.skin = 0.4
        self.p0 = self.system.part.add(
            pos=[0, 0, 0], type=0, fix=3 * [True], q=self.q1)
        self.p1 = self.system.part.add(
            pos=[2, 0, 0], type=0, fix=3 * [True], q=self.q2)

        p3m = espressomd.electrostatics.P3M(
            prefactor=COULOMB_PREFACTOR, accuracy=1e-6, mesh=3 * [52], cao=4)
        self.system.electrostatics.solver = p3m

        self.system.non_bonded_inter[0, 0].thole.set_params(
            scaling_coeff=self.thole_s, q1q2=self.q1 * self.q2)

    def tearDown(self):
        self.system.electrostatics.solver = None
        self.system.part.clear()

    def test(self):
        ns = 100
        for i in range(1, ns):
            x = 20.0 * i / ns
            self.p1.pos = [x, 0, 0]
            self.system.integrator.run(0)

            sd = x * self.thole_s
            # Force is exact
            F_calc = COULOMB_PREFACTOR * self.q1 * self.q2 / x**2 * \
                0.5 * (2.0 - (np.exp(-sd) * (sd * (sd + 2.0) + 2.0)))
            # Energy is slightly off due to self-energy.
            # Error is approximated with erfc for given system parameters
            E_calc = COULOMB_PREFACTOR * self.q1 * self.q2 / x * \
                (1.0 - np.exp(-sd) * (1.0 + sd / 2.0)) - \
                0.250088 * math.erfc(0.741426 * x)

            E = self.system.analysis.energy()
            self.assertAlmostEqual(self.p1.f[0], F_calc, delta=1e-3)
            self.assertAlmostEqual(E["total"], E_calc, delta=1.2e-2)


if __name__ == "__main__":
    ut.main()
