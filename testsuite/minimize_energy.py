from __future__ import print_function
import sys
import unittest as ut
import numpy as np
import espressomd
from espressomd.interactions import *

@ut.skipIf(not espressomd.has_features("LENNARD_JONES"),"Skipped because LENNARD_JONES turned off.")
class test_minimize_energy(ut.TestCase):
    system = espressomd.System()
    box_l = 10.0
    density = 0.6
    vol = box_l * box_l *box_l
    n_part = int( vol * density )

    lj_eps = 1.0
    lj_sig = 1.0
    lj_cut = 1.12246


    def runTest(self):
        self.system.box_l = [self.box_l, self.box_l, self.box_l]
        self.system.cell_system.skin = 0.4
        self.system.time_step = 0.01
        self.system.non_bonded_inter[0, 0].lennard_jones.set_params(
            epsilon=self.lj_eps, sigma=self.lj_sig,
            cutoff=self.lj_cut, shift="auto")

        for i in range(self.n_part):
            self.system.part.add(id=i, pos=np.random.random(3) * self.system.box_l)

        self.system.integrator.set_steepest_descent(f_max=0.0, gamma=0.1, max_displacement=0.001)


        self.system.integrator.run(10000)
        self.system.integrator.run(10000)

        energy = self.system.analysis.energy()

        self.assertEqual(energy["total"], 0)

if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
