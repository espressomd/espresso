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

import espressomd
from espressomd import lb, lbboundaries, shapes, has_features
import unittest as ut
import numpy as np
import sys

# Define the LB Parameters
TIME_STEP = 0.1
AGRID = 1.0 
KVISC = 7 
DENS = 1 
BOX_SIZE = 6 * AGRID
F = 1. / BOX_SIZE**3

LB_PARAMS = {'agrid': AGRID,
             'dens': DENS,
             'visc': KVISC,
             'tau': TIME_STEP,
             'ext_force_density': [0, F, 0]}


class Momentum(object):
    lbf = None
    system = espressomd.System(box_l=[BOX_SIZE] * 3)
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.01

    def test(self):
        print(self.system.cell_system.get_state())
        self.system.actors.clear()
        self.system.part.clear()
        self.system.actors.add(self.lbf)
        self.system.thermostat.set_lb(LB_fluid=self.lbf, gamma=1, seed=1)

        applied_force = self.system.volume() * np.array(
            LB_PARAMS['ext_force_density'])
        p = self.system.part.add(pos=(0, 0, 0), ext_force=-applied_force)

        # Reach steady state
        self.system.integrator.run(3000)
        v_final = np.copy(p.v)

        for i in range(30):
            self.system.integrator.run(1000)
            np.testing.assert_allclose(
                self.system.analysis.linear_momentum(),
              [0, 0, 0], atol=1E-8)
#            np.testing.assert_allclose(
#                p.v * p.mass, -
#                    np.array(
#                        self.system.analysis.linear_momentum(
#                            include_particles=False))
#                            +np.array((0,F *AGRID/ TIME_STEP /2 ,0)),
#                atol=np.linalg.norm(applied_force) * TIME_STEP * 0.55)

            # Check that particle velocity is stationary
            # up to the acceleration of 1/2 time step
            np.testing.assert_allclose(np.copy(p.v), v_final, 
                                       atol=np.linalg.norm(applied_force) / p.mass * TIME_STEP * 0.55)
           

@ut.skipIf(not espressomd.gpu_available() or not espressomd.has_features(
    ['EXTERNAL_FORCES']), "Skipping test due to missing features.")
class LBGPUMomentum(ut.TestCase, Momentum):

    def setUp(self):
        self.lbf = espressomd.lb.LBFluidGPU(**LB_PARAMS)


@ut.skipIf(not espressomd.has_features(
    ['EXTERNAL_FORCES']), "Skipping test due to missing features.")
class LBCPUMomentum(ut.TestCase, Momentum):

    def setUp(self):
        self.lbf = espressomd.lb.LBFluid(**LB_PARAMS)


@ut.skipIf(not espressomd.has_features(
    ['LB_WALBERLA', 'EXTERNAL_FORCES']), "Skipping test due to missing features.")
class LBWalberlaMomentum(ut.TestCase, Momentum):

    def setUp(self):
        self.lbf = espressomd.lb.LBFluidWalberla(**LB_PARAMS)


if __name__ == "__main__":
    ut.main()
