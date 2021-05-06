# Copyright (C) 2010-2020 The ESPResSo project
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

import espressomd.lb
from thermostats_common import ThermostatsCommon

"""
Check the lattice-Boltzmann thermostat with respect to the particle velocity
distribution.


"""

KT = 0.9 
AGRID = 0.8
node_volume = AGRID**3
VISC = 6 
DENS = 1.7
TIME_STEP = 0.005
GAMMA = 2
LB_PARAMS = {'agrid': AGRID,
             'dens': DENS,
             'visc': VISC,
             'tau': TIME_STEP,
             'kT': KT,
             'seed': 123}


class LBThermostatCommon(ThermostatsCommon):

    """Base class of the test that holds the test logic."""
    lbf = None
    system = espressomd.System(box_l=[AGRID * 12] * 3)
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.4 * AGRID

    def test_fluid(self):
        self.prepare()
        self.system.integrator.run(100)
        fluid_temps = []
        for _ in range(100):
            fluid_temps.append(
                np.average([n.density * n.velocity**2 for n in self.lbf.nodes()]) * node_volume)
            self.system.integrator.run(3)

        fluid_temp = np.average(fluid_temps)
        self.assertAlmostEqual(fluid_temp, KT, delta=0.05)

    def prepare(self):
        self.system.actors.clear()
        self.system.actors.add(self.lbf)
        self.system.thermostat.set_lb(LB_fluid=self.lbf, seed=5, gamma=GAMMA)

    def test_with_particles(self):
        self.prepare()
        self.system.part.add(
            pos=np.random.random((100, 3)) * self.system.box_l)
        self.system.integrator.run(100)
        N = len(self.system.part)
        loops = 500
        v_particles = np.zeros((loops, N, 3))
        fluid_temps = []

        for i in range(loops):
            self.system.integrator.run(3)
            if i % 10 == 0:
                fluid_temps.append(
                    np.average([n.density * n.velocity**2 for n in self.lbf.nodes()]) * node_volume)
            v_particles[i] = self.system.part[:].v
        fluid_temp = np.average(fluid_temps)

        np.testing.assert_allclose(np.average(v_particles), 0, atol=0.02)
        np.testing.assert_allclose(np.var(v_particles), KT, rtol=0.02)

        minmax = 3
        n_bins = 7
        error_tol = 0.015  
        self.check_velocity_distribution(
            v_particles.reshape((-1, 3)), minmax, n_bins, error_tol, KT)

        np.testing.assert_allclose(fluid_temp, KT, rtol=0.01)


@utx.skipIfMissingFeatures("LB_WALBERLA")
class LBWalberlaThermostat(ut.TestCase, LBThermostatCommon):

    """Test for the CPU implementation of the LB."""

    def setUp(self):
        self.lbf = espressomd.lb.LBFluidWalberla(**LB_PARAMS)


if __name__ == '__main__':
    ut.main()
