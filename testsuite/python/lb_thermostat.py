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
import thermostats_common
import numpy as np

import espressomd.lb

"""
Check the lattice-Boltzmann thermostat with respect to the particle velocity
distribution.
"""

KT = 0.25
AGRID = 2.5
VISC = 2.7
DENS = 1.7
TIME_STEP = 0.05
LB_PARAMS = {'agrid': AGRID,
             'dens': DENS,
             'visc': VISC,
             'tau': TIME_STEP,
             'kT': KT,
             'seed': 123}


class LBThermostatCommon(thermostats_common.ThermostatsCommon):

    """Check the LB thermostat."""

    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.4 * AGRID
    np.random.seed(42)

    def setUp(self):
        self.lbf = self.lb_class(**LB_PARAMS)
        self.system.actors.add(self.lbf)

    def tearDown(self):
        self.system.actors.clear()
        self.system.thermostat.turn_off()

    def test_velocity_distribution(self):
        self.system.part.add(
            pos=np.random.random((100, 3)) * self.system.box_l)
        self.system.thermostat.set_lb(LB_fluid=self.lbf, seed=5, gamma=5.0)
        self.system.integrator.run(20)
        N = len(self.system.part)
        loops = 250
        v_stored = np.zeros((loops, N, 3))
        for i in range(loops):
            self.system.integrator.run(3)
            v_stored[i] = self.system.part[:].v
        minmax = 5
        n_bins = 7
        error_tol = 0.01
        self.check_velocity_distribution(
            v_stored.reshape((-1, 3)), minmax, n_bins, error_tol, KT)


class LBCPUThermostat(LBThermostatCommon, ut.TestCase):

    """Test for the CPU implementation of the LB."""

    lb_class = espressomd.lb.LBFluid


@utx.skipIfMissingGPU()
class LBGPUThermostat(LBThermostatCommon, ut.TestCase):

    """Test for the GPU implementation of the LB."""

    lb_class = espressomd.lb.LBFluidGPU


if __name__ == '__main__':
    ut.main()
