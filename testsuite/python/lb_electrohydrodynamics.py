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
import espressomd.lb


class LBEHTest(object):
    from espressomd import lb
    s = espressomd.System(box_l=[6.0, 6.0, 6.0])

    def setUp(self):
        self.params = {'time_step': 0.01,
                       'tau': 0.02,
                       'agrid': 0.5,
                       'dens': 0.85,
                       'viscosity': 30.0,
                       'friction': 3.0,
                       'temp': 0.0,
                       'skin': 0.2,
                       'muE': [0.1, 0.2, 0.3]}

        self.s.periodicity = [1, 1, 1]
        self.s.time_step = self.params['time_step']
        self.s.cell_system.skin = self.params['skin']

        self.lbf = self.LBClass(
            visc=self.params['viscosity'],
            dens=self.params['dens'],
            agrid=self.params['agrid'],
            tau=self.s.time_step,
            kT=self.params['temp']
        )

        self.s.actors.add(self.lbf)
        self.s.thermostat.set_lb(
            LB_fluid=self.lbf,
            gamma=self.params['friction'])

    def tearDown(self):
        self.s.actors.clear()

    def test(self):
        s = self.s

        p = s.part.add(pos=0.5 * self.s.box_l, mu_E=self.params['muE'])

        mu_E = np.array(self.params['muE'])
        # Terminal velocity is mu_E minus the momentum the fluid
        # got by accelerating the particle in the beginning.
        v_term = (1. - 1. / (s.volume() * self.params['dens'])) * mu_E

        s.integrator.run(steps=500)

        np.testing.assert_allclose(v_term, np.copy(p.v), atol=5e-5)


@utx.skipIfMissingFeatures(["LB_WALBERLA", "LB_ELECTROHYDRODYNAMICS"])
class LBEHWalberla(LBEHTest, ut.TestCase):

    def setUp(self):
        self.LBClass = espressomd.lb.LBFluidWalberla
        LBEHTest.setUp(self)


if __name__ == "__main__":
    ut.main()
