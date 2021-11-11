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

import espressomd
import unittest as ut
import unittest_decorators as utx
import numpy as np

# Define the LB parameters
TIME_STEP = 0.1
AGRID = 1.0
KVISC = 5
DENS = 1
BOX_SIZE = 6 * AGRID
F = 1. / BOX_SIZE**3
GAMMA = 15

LB_PARAMS = {'agrid': AGRID,
             'dens': DENS,
             'visc': KVISC,
             'tau': TIME_STEP,
             'ext_force_density': [0, F, 0]}


class Momentum(object):
    """
    Tests momentum conservation for an LB coupled to a particle, where opposing
    forces are applied to LB and particle. The test should uncover issues
    with boundary and ghost layer handling.

    """
    system = espressomd.System(box_l=[BOX_SIZE] * 3)
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.01

    def setUp(self):
        self.lbf = self.lb_class(**LB_PARAMS)

    def tearDown(self):
        self.system.actors.clear()
        self.system.thermostat.turn_off()
        self.system.part.clear()

    def test(self):
        self.system.actors.add(self.lbf)
        self.system.thermostat.set_lb(LB_fluid=self.lbf, gamma=GAMMA, seed=1)
        applied_force = self.system.volume() * np.array(
            LB_PARAMS['ext_force_density'])
        p = self.system.part.add(
            pos=(0, 0, 0), ext_force=-applied_force, v=[.1, .2, .3])

        # Reach steady state
        self.system.integrator.run(500)
        v_final = np.copy(p.v)
        momentum = self.system.analysis.linear_momentum()

        for _ in range(10):
            self.system.integrator.run(50)
            # check that momentum stays constant
            np.testing.assert_allclose(
                self.system.analysis.linear_momentum(), momentum, atol=2E-4)

            # Check that particle velocity is stationary
            # up to the acceleration of 1/2 time step
            np.testing.assert_allclose(np.copy(p.v), v_final, atol=2.2E-3)

        # Make sure, the particle has crossed the periodic boundaries
        self.assertGreater(
            np.amax(
                np.abs(v_final) *
                self.system.time),
            BOX_SIZE)


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(['EXTERNAL_FORCES'])
class LBGPUMomentum(Momentum, ut.TestCase):

    lb_class = espressomd.lb.LBFluidGPU


@utx.skipIfMissingFeatures(['EXTERNAL_FORCES'])
class LBCPUMomentum(Momentum, ut.TestCase):

    lb_class = espressomd.lb.LBFluid


if __name__ == "__main__":
    ut.main()
