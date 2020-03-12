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
import numpy as np
import espressomd


@utx.skipIfMissingFeatures(["ENGINE"])
class SwimmerTest(ut.TestCase):

    def test(self):
        boxl = 12
        sampsteps = 2000
        tstep = 0.01

        v_swim = 0.3
        f_swim = 0.1
        temp = 0.0
        gamma = 1.0

        pos_0 = np.array([boxl / 2., boxl / 2., 1. * boxl / 3.])
        pos_1 = np.array([boxl / 2., boxl / 2., 2. * boxl / 3.])

        def z_f(t, z0):
            return f_swim / gamma * \
                (-1. / gamma + t + (1. / gamma) * np.exp(-gamma * t)) + z0

        def z_v(t, z0):
            return v_swim * (-1. / gamma + t + (1. / gamma) *
                             np.exp(-gamma * t)) + z0

        S = espressomd.System(box_l=[1.0, 1.0, 1.0])

        S.box_l = [boxl, boxl, boxl]
        S.cell_system.skin = 0.1
        S.time_step = tstep

        S.part.add(id=0, pos=pos_0, swimming={"v_swim": v_swim})
        S.part.add(id=1, pos=pos_1, swimming={"f_swim": f_swim})
        S.part[:].rotation = 1, 1, 1

        S.thermostat.set_langevin(kT=temp, gamma=gamma, seed=42)

        S.integrator.run(sampsteps)

        pos_0[2] = z_v(S.time, pos_0[2])
        pos_1[2] = z_f(S.time, pos_1[2])

        delta_pos_0 = np.linalg.norm(S.part[0].pos - pos_0)
        delta_pos_1 = np.linalg.norm(S.part[1].pos - pos_1)

        self.assertLess(1.4e-3, delta_pos_0)
        self.assertLess(delta_pos_0, 1.6e-3)
        self.assertLess(4.9e-4, delta_pos_1)
        self.assertLess(delta_pos_1, 5.1e-4)


if __name__ == '__main__':
    ut.main()
