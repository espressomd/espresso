from __future__ import print_function
import unittest as ut
import espressomd
import espressomd.analyze
import espressomd.lb
import numpy as np

@ut.skipIf((not espressomd.has_features(["LB_GPU"])) or
           espressomd.has_features(["SHANCHEN"]), "Features not available, skipping test!")
class RemoveTotalMomentumTest(ut.TestCase):
    def test(self):
        dt = 0.01
        skin = 0.1
        agrid = 1.0
        fric = 20.0
        visc = 1.0
        dens = 1.0

        s = espressomd.System()
        s.box_l = [10, 10, 10]
        s.time_step = dt
        s.cell_system.skin = skin

        for i in range(100):
            r = s.box_l * np.random.random(3)
            v = [0., 0., 1.]
            s.part.add(pos=r, v=v)

        lbf = espressomd.lb.LBFluid_GPU(
            agrid=agrid, fric=fric, dens=dens, visc=visc, tau=dt)

        s.actors.add(lbf)

        s.integrator.run(300)

        lbf.remove_total_momentum()

        p = np.array(s.analysis.analyze_linear_momentum())

        self.assertTrue(np.all(np.abs(p) < 1e-3))


if __name__ == "__main__":
    #print("Features: ", espressomd.features())
    ut.main()
