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
        dens = 12.0

        s = espressomd.System(box_l=[1.0, 1.0, 1.0])
        s.seed = s.cell_system.get_state()['n_nodes'] * [1234]
        s.box_l = [10, 10, 10]
        s.time_step = dt
        s.cell_system.skin = skin

        for i in range(100):
            r = s.box_l * np.random.random(3)
            v = [1., 1., 1.] * np.random.random(3)
            # Make sure that id gaps work correctly
            s.part.add(id=2*i, pos=r, v=v)

        if espressomd.has_features(["MASS"]):
            # Avoid masses too small for the time step
            s.part[:].mass = 2. * (0.1 + np.random.random(100))

        lbf = espressomd.lb.LBFluidGPU(
            agrid=agrid, fric=fric, dens=dens, visc=visc, tau=dt)

        s.actors.add(lbf)

        s.integrator.run(300)

        lbf.remove_total_momentum()

        p = np.array(s.analysis.analyze_linear_momentum())

        self.assertAlmostEqual(np.max(p), 0., places=3)
        self.assertAlmostEqual(np.min(p), 0., places=3)

if __name__ == "__main__":
    #print("Features: ", espressomd.features())
    ut.main()
