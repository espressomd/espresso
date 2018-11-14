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
import unittest as ut
import numpy as np

import espressomd
import espressomd.lb

"""
Check the Lattice Boltzmann mass conservation.


"""

KT = 2.25
AGRID = .5
VISC = .7
DENS = 1.7
TIME_STEP = 0.01
LB_PARAMS = {'agrid': AGRID,
             'dens': DENS,
             'visc': VISC,
             'fric': 2.0,
             'tau': TIME_STEP}


class LBThermostatCommon(object):

    """Base class of the test that holds the test logic."""
    lbf = None
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.4 * AGRID

    def prepare(self):
        self.system.set_random_state_PRNG()
        self.system.actors.clear()
        self.system.actors.add(self.lbf)
        self.system.thermostat.set_lb(kT=KT)

    def test_mass_conservation(self):
        self.prepare()
        self.system.integrator.run(1000)
        nodes=[]
        for i in range(int(self.system.box_l[0]/LB_PARAMS['agrid'])):
          for j in range(int(self.system.box_l[1]/LB_PARAMS['agrid'])):
            for k in range(int(self.system.box_l[2]/LB_PARAMS['agrid'])):
              nodes.append(self.lbf[i,j,k])
          
        result = np.zeros(10)
        for i in range(10):
          self.system.integrator.run(10)
          v=[]
          for n in nodes:
            v.append(n.density[0])
          result[i] = np.mean(v)
        np.testing.assert_allclose(result - DENS, np.zeros_like(result), atol=1e-7)
        

@ut.skipIf(not espressomd.has_features(
    ['LB']) or espressomd.has_features("SHANCHEN"), "Skipping test due to missing features.")
class LBCPUThermostat(ut.TestCase, LBThermostatCommon):

    """Test for the CPU implementation of the LB."""

    def setUp(self):
        self.lbf = espressomd.lb.LBFluid(**LB_PARAMS)


@ut.skipIf(not espressomd.has_features(
    ['LB_GPU']) or espressomd.has_features("SHANCHEN"), "Skipping test due to missing features.")
class LBGPUThermostat(ut.TestCase, LBThermostatCommon):

    """Test for the GPU implementation of the LB."""

    def setUp(self):
        self.lbf = espressomd.lb.LBFluidGPU(**LB_PARAMS)


if __name__ == '__main__':
    ut.main()
