#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
#
from __future__ import print_function
import unittest as ut
import espressomd
import numpy as np
from espressomd import electrostatics
from tests_common import *


@ut.skipIf(not espressomd.has_features(["ELECTROSTATICS", "EXTERNAL_FORCES"]),
           "Features not available, skipping test!")
class ElectrostaticInteractionsTests(ut.TestCase):
    # Handle to espresso system
    system = espressomd.System()
    system.seed = system.cell_system.get_state()['n_nodes'] * [1234]
    p3m_energy = -0.501062398379
    p3m_force = 2.48921612e-01

    def setUp(self):
        self.system.box_l = [20, 20, 20]
        self.system.time_step = 0.01

        if not self.system.part.exists(0):
            self.system.part.add(id=0, pos=(1.0, 2.0, 2.0), q=1, fix=[1,1,1])
        if not self.system.part.exists(1):
            self.system.part.add(id=1, pos=(3.0, 2.0, 2.0), q=-1, fix=[1,1,1])
        print("ut.TestCase setUp")

    @ut.skipIf(not espressomd.has_features(["P3M"]),
               "Features not available, skipping test!")
    def test_p3m(self):

        test_P3M = generate_test_for_class(
            self.system,
            electrostatics.P3M,
            dict(
                 accuracy=9.910945054074526e-08,
                 mesh=[22,22,22],
                 cao=7,
                 r_cut=8.906249999999998,
                 alpha=0.387611049779351,
                 tune=False))
        p3m = espressomd.electrostatics.P3M(bjerrum_length=1.0,
                                            accuracy=9.910945054074526e-08,
                                            mesh=[22,22,22],
                                            cao=7,
                                            r_cut=8.906249999999998,
                                            alpha=0.387611049779351,
                                            tune=False)
        self.system.actors.add(p3m)
        
        self.assertAlmostEqual(self.system.analysis.energy()['coulomb'],
                               self.p3m_energy)

        # need to update forces
        self.system.integrator.run(steps=1, recalc_forces=True)
        self.assertTrue( np.allclose(self.system.part[0].f,
                                     [self.p3m_force, 0, 0]))
        self.assertTrue( np.allclose(self.system.part[1].f,
                                     [-self.p3m_force, 0, 0]))
        

    @ut.skipIf(not espressomd.has_features(["DEBYE_HUECKEL"]),
               "Features not available, skipping test!")
    def test_dh(self):
        
        test_CDH = generate_test_for_class(
            elf.system,
            electrostatics.CDH,
            dict(
                bjerrum_length=1.0,
                kappa=2.3,
                r_cut=2,
                r0=1,
                r1=1.9,
                eps_int=0.8,
                eps_ext=1,
                alpha=2))



        #test_DH = generate_test_for_class(
            #self.system,
            #electrostatics.P3M,
            #dict(
                #bjerrum_length=1.0,
                #kappa=2.0,
                #r_cut=2))
            
        #dh = espressomd.electrostatics.CDH(bjerrum_length=1.0,
                                           #kappa=2.0,
                                           #r_cut=2)
        #self.system.actors.add(dh)
        
        #self.assertAlmostEqual(self.system.analysis.energy()['coulomb'],
                               #self.p3m_energy)

        ## need to update forces
        #self.system.integrator.run(steps=1, recalc_forces=True)
        #self.assertTrue( np.allclose(self.system.part[0].f,
                                     #[self.p3m_force, 0, 0]))
        #self.assertTrue( np.allclose(self.system.part[1].f,
                                     #[-self.p3m_force, 0, 0]))




if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
