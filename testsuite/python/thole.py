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
from __future__ import print_function
import unittest as ut
import espressomd
import numpy as np
import espressomd.interactions

COULOMB_PREFACTOR = 1.3


@ut.skipIf(not espressomd.has_features(["THOLE", "EXTERNAL_FORCES", "P3M"]),
           "Features not available, skipping test!")
class TestThole(ut.TestCase):

    """
    This testcase takes a large box to minimize periodic effects and tests the
    thole damping nonbonded interaction forces agains the analytical result
    """

    box_l = 500.0
    
    thole_s = 1.0
    q1 = 1.0
    q2 = -1.0

    system = espressomd.System(box_l=[box_l] * 3)

    def setUp(self):
        from espressomd.electrostatics import P3M

        self.system.time_step = 0.01
        self.system.cell_system.skin = 0.4


        self.system.part.add(pos=[0, 0, 0], type=0, fix=[1, 1, 1], q=self.q1)
        self.system.part.add(pos=[2, 0, 0], type=0, fix=[1, 1, 1], q=self.q2)

        p3m = P3M(
            prefactor=COULOMB_PREFACTOR,
            accuracy=1e-6,
            mesh=[52,
                  52,
                  52],
            cao=4)
        self.system.actors.add(p3m)

        self.system.non_bonded_inter[0, 0].thole.set_params(
            scaling_coeff=self.thole_s, q1q2=self.q1 * self.q2)

    def test(self):
        res_dForce = []
        res_dEnergy = []
        Es=[]
        ns = 100
        for i in range(1, ns):
            x = 8.0 * i / ns
            self.system.part[1].pos = [x, 0, 0]
            self.system.integrator.run(0)
            res_dForce.append(self.system.part[1].f[0] - (COULOMB_PREFACTOR * self.q1*self.q2 / x**2 * 0.5 *
                                                   (2.0 - (np.exp(-x) * (x * (x + 2.0) + 2.0)))))
            res_dEnergy.append(self.system.analysis.energy()["total"] - (COULOMB_PREFACTOR * self.q1*self.q2 / x *
                                                   (1.0 - np.exp(-self.thole_s * x) * (1.0 + self.thole_s * x /2.0))))

            Es.append([self.system.analysis.energy()["total"], COULOMB_PREFACTOR * self.q1*self.q2 / x *
                                                   (1.0 - np.exp(-self.thole_s * x) * (1.0 + self.thole_s * x /2.0))])

        np.savetxt('thole.data',Es)
        for f in res_dForce:
            self.assertLess(
                abs(f), 1e-3, msg="Deviation of thole interaction force (damped coulomb) from analytical result too large")
        
        for e in res_dEnergy:
            self.assertLess(
                abs(e), 1e-3, msg="Deviation of thole interaction energy (damped coulomb) from analytical result too large")


if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
