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
TIME_STEP = 0.008
AGRID = .4 
GRID_SIZE = 6 
KVISC = 4
DENS = 2.3
F = 5.5 / GRID_SIZE**3 
GAMMA = 1


LB_PARAMS = {'agrid': AGRID,
             'dens': DENS,
             'visc': KVISC,
             'tau': TIME_STEP,
             'ext_force_density': [-.7 * F, .9 * F, .8 * F]}


class Momentum(object):
    """
    Tests momentum conservation for an LB coupled to a particle, where opposing
    forces are applied to LB and particle. The test should uncover issues
    with boundary and ghost layer handling.

    """
    lbf = None
    system = espressomd.System(box_l=[GRID_SIZE * AGRID] * 3)
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.01

    def test(self):
        self.system.actors.clear()
        self.system.part.clear()
        self.system.actors.add(self.lbf)
        self.system.thermostat.set_lb(LB_fluid=self.lbf, gamma=GAMMA, seed=1)
        np.testing.assert_allclose(
            self.lbf.ext_force_density,
            LB_PARAMS["ext_force_density"])

        # Initial momentum before integration = 0
        np.testing.assert_allclose(
            self.system.analysis.linear_momentum(), [0., 0., 0.], atol=1E-12)

        ext_fluid_force = self.system.volume() * np.array(
            LB_PARAMS['ext_force_density'])

        p = self.system.part.add(
            pos=self.system.box_l / 2, ext_force=-ext_fluid_force, v=[.2, .4, .6])
        initial_momentum = np.array(self.system.analysis.linear_momentum())
        np.testing.assert_allclose(initial_momentum, np.copy(p.v) * p.mass)
        while True: 
            self.system.integrator.run(500)

            measured_momentum = self.system.analysis.linear_momentum()
            coupling_force = -(p.f - p.ext_force)
            compensation = -TIME_STEP / 2 * coupling_force

            np.testing.assert_allclose(measured_momentum + compensation, 
                                       initial_momentum, atol=1E-4)
            if np.linalg.norm(p.f) < 0.01 \
               and np.all(np.abs(p.pos) > 10.1 * self.system.box_l):
                break

        # Make sure, the particle has crossed the periodic boundaries
        self.assertGreater(
            max(
                np.abs(p.v) *
                self.system.time),
            self.system.box_l[0])


@utx.skipIfMissingFeatures(['LB_WALBERLA', 'EXTERNAL_FORCES'])
class LBWalberlaMomentum(ut.TestCase, Momentum):

    def setUp(self):
        self.lbf = espressomd.lb.LBFluidWalberla(**LB_PARAMS)


if __name__ == "__main__":
    ut.main()
