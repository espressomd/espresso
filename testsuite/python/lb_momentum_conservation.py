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

# Define the LB Parameters
TIME_STEP = 0.005
AGRID = 1.0
KVISC = 5
DENS = 1
BOX_SIZE = 6 * AGRID
F = 0.05
GAMMA = 5

LB_PARAMS = {'agrid': AGRID,
             'dens': DENS,
             'visc': KVISC,
             'tau': TIME_STEP,
             'ext_force_density': [-.5 * F, F, .8 *F]}


class Momentum(object):
    """
    Tests momentum conservation for an LB coupled to a particle, where opposing
    forces are applied to LB and particle. The test should uncover issues
    with boundary and ghost layer handling.

    """
    lbf = None
    system = espressomd.System(box_l=[BOX_SIZE] * 3)
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.01

    def test(self):
        self.system.actors.clear()
        self.system.part.clear()
        self.system.actors.add(self.lbf)
        self.system.thermostat.set_lb(LB_fluid=self.lbf, gamma=GAMMA, seed=1)
        np.testing.assert_allclose(self.lbf.ext_force_density, LB_PARAMS["ext_force_density"])
        
        # Initial momentum before integration = 0
        np.testing.assert_allclose(
          self.system.analysis.linear_momentum(), [0., 0., 0.])
        
        applied_force = self.system.volume() * np.array(
            LB_PARAMS['ext_force_density'])

        p = self.system.part.add(
            pos=(0, 0, 0), ext_force=-applied_force, v=[.1, .2, .3])
        initial_momentum = self.system.analysis.linear_momentum()
        np.testing.assert_allclose(initial_momentum, np.copy(p.v) * p.mass)
        
        while True: 
            self.system.integrator.run(100)
            measured_momentum = self.system.analysis.linear_momentum()
            f_half_compensation = self.lbf.ext_force_density * TIME_STEP *DENS *self.system.volume()/2
            stokes_force_compensation = p.f * TIME_STEP  # fluid force is opposed to particle force
            np.testing.assert_allclose(measured_momentum + f_half_compensation +stokes_force_compensation, initial_momentum, atol=1.4E-4)
            if np.linalg.norm(stokes_force_compensation) < 1E-4: break


        # Make sure, the particle has crossed the periodic boundaries
        self.assertGreater(
            max(
                np.abs(p.v) *
                self.system.time),
            BOX_SIZE)


@utx.skipIfMissingFeatures(['LB_WALBERLA', 'EXTERNAL_FORCES'])
class LBWalberlaMomentum(ut.TestCase, Momentum):

    def setUp(self):
        self.lbf = espressomd.lb.LBFluidWalberla(**LB_PARAMS)


if __name__ == "__main__":
    ut.main()
