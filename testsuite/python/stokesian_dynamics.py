import espressomd
from espressomd import constraints
from espressomd.thermostat import flags
import numpy as np
import matplotlib.pyplot as plt
import unittest as ut
import unittest_decorators as utx
from tests_common import abspath

@utx.skipIfMissingFeatures(["STOKESIAN_DYNAMICS"])
class StokesianDynamicsTest(ut.TestCase):
    system = espressomd.System(box_l=[1.0,1.0,1.0])
    data = np.loadtxt(abspath('data/dancing.txt'))

    def setUp(self):
        self.system.part.clear()
        self.system.box_l = [10] * 3

        self.system.time_step = 1.0
        self.system.cell_system.skin = 0.4

        self.system.part.add(pos=[-5,0,0],rotation=[1,1,1])
        self.system.part.add(pos=[ 0,0,0],rotation=[1,1,1])
        self.system.part.add(pos=[ 7,0,0],rotation=[1,1,1])

        self.system.thermostat.set_sd(viscosity=1.0,
                                      device="cpu",
                                      radii={ 0: 1.0 },
                                      flags=flags.SELF_MOBILITY | flags.PAIR_MOBILITY | flags.FTS)
        self.system.integrator.set_sd()

        gravity = constraints.Gravity(g=[0,-1,0])
        self.system.constraints.add(gravity)

    def test(self):
        intsteps = 8000
        self.pos = np.empty([intsteps+1,3*len(self.system.part)])

        self.pos[0,:] = self.system.part[:].pos.flatten()
        for i in range(intsteps):
            self.system.integrator.run(1)
            for n,p in enumerate(self.system.part):
                self.pos[i+1,3*n:3*n+3] = p.pos

        for i in range(self.data.shape[0]):
            for n in range(3):
                x_ref = self.data[i,2*n]
                y_ref = self.data[i,2*n+1]
                x = self.pos[:,3*n]
                y = self.pos[:,3*n+1]

                if y_ref < -555:
                    continue

                idx = np.abs(y - y_ref).argmin()
                self.assertTrue(np.abs(y_ref - y[idx]) < 0.5)

if __name__ == '__main__':
    ut.main()
