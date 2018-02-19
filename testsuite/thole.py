from __future__ import print_function
import unittest as ut
import espressomd
import numpy as np
import espressomd.interactions

@ut.skipIf(not espressomd.has_features(["THOLE", "EXTERNAL_FORCES"]),
           "Features not available, skipping test!")

class TestThole(ut.TestCase):
    """
    This testcase takes a large box to minimize periodic effects and tests the
    thole damping nonbonded interaction forces agains the analytical result
    """

    box_l = 500.0
    system = espressomd.System(box_l = [box_l]*3)

    def setUp(self):
        from espressomd.electrostatics import P3M

        self.system.time_step = 0.01
        self.system.cell_system.skin = 0.4

        q1 = 1.0
        q2 = -1.0

        self.system.part.add(pos=[0, 0, 0], type=0, fix=[1, 1, 1], q=q1)
        self.system.part.add(pos=[2, 0, 0], type=0, fix=[1, 1, 1], q=q2)

        p3m = P3M(prefactor=1.0, accuracy=1e-6, mesh=[52, 52, 52], cao=4)
        self.system.actors.add(p3m)

        self.system.non_bonded_inter[0, 0].thole.set_params(
            scaling_coeff=1.0, q1q2=q1 * q2)

    def test(self):
        res = []
        ns = 100
        for i in range(1, ns):
            x = 8.0 * i / ns
            self.system.part[1].pos = [x, 0, 0]
            self.system.integrator.run(0)
            res.append(self.system.part[1].f[0] - (-1 / x**2 * 0.5 *
                                         (2.0 - (np.exp(-x) * (x * (x + 2.0) + 2.0)))))

        for f in res:
            self.assertLess(abs(f),1e-3, msg = "Deviation of thole interaction (damped coulomb) from analytical result too large")


if __name__ == "__main__":
    print("Features: ", espressomd.features())
    ut.main()
