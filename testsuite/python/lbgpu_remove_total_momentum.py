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
import espressomd.analyze
import espressomd.lb
import numpy as np


@ut.skipIf((not espressomd.gpu_available() or not espressomd.has_features(["CUDA"])), "Features or gpu not available, skipping test!")
class RemoveTotalMomentumTest(ut.TestCase):

    def test(self):
        dt = 0.01
        skin = 0.1
        agrid = 1.0
        fric = 20.0
        visc = 1.0
        dens = 12.0

        s = espressomd.System(box_l=[1.0, 1.0, 1.0])
        s.seed = s.cell_system.get_state()['n_nodes'] * [1234]
        s.box_l = [10, 10, 10]
        s.time_step = dt
        s.cell_system.skin = skin

        for i in range(100):
            r = s.box_l * np.random.random(3)
            v = [1., 1., 1.] * np.random.random(3)
            # Make sure that id gaps work correctly
            s.part.add(id=2 * i, pos=r, v=v)

        if espressomd.has_features(["MASS"]):
            # Avoid masses too small for the time step
            s.part[:].mass = 2. * (0.1 + np.random.random(100))

        lbf = espressomd.lb.LBFluidGPU(
            agrid=agrid, dens=dens, visc=visc, tau=dt)
        s.thermostat.set_lb(LB_fluid=lbf, gamma=fric)
        s.actors.add(lbf)

        s.integrator.run(300)

        lbf.remove_total_momentum()

        p = np.array(s.analysis.analyze_linear_momentum())

        self.assertAlmostEqual(np.max(p), 0., places=3)
        self.assertAlmostEqual(np.min(p), 0., places=3)

if __name__ == "__main__":
    ut.main()
