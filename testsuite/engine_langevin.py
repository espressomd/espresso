from __future__ import print_function
import unittest as ut
import numpy as np
import espressomd


@ut.skipIf(not espressomd.has_features(["ENGINE"]),
           "Features not available, skipping test!")
class SwimmerTest(ut.TestCase):
    def test(self):
        boxl = 12
        sampsteps = 2000
        tstep = 0.01

        v_swim = 0.3
        f_swim = 0.1
        temp = 0.0
        gamma = 1.0

        pos_0 = np.array([boxl / 2., boxl / 2., 1. * boxl / 3.])
        pos_1 = np.array([boxl / 2., boxl / 2., 2. * boxl / 3.])

        def z_f(t, z0):
            return f_swim / gamma * \
                (-1. / gamma + t + (1. / gamma) * np.exp(-gamma * t)) + z0

        def z_v(t, z0):
            return v_swim * (-1. / gamma + t + (1. / gamma) *
                             np.exp(-gamma * t)) + z0

        S = espressomd.System(box_l=[1.0, 1.0, 1.0])
        S.seed  = S.cell_system.get_state()['n_nodes'] * [1234]

        S.box_l = [boxl, boxl, boxl]
        S.cell_system.skin = 0.1
        S.time_step = tstep

        S.part.add(id=0, pos=pos_0, swimming={"v_swim": v_swim})
        S.part.add(id=1, pos=pos_1, swimming={"f_swim": f_swim})
        S.part[:].rotation=1,1,1

        S.thermostat.set_langevin(kT=temp, gamma=gamma)

        S.integrator.run(sampsteps)

        pos_0[2] = z_v(S.time, pos_0[2])
        pos_1[2] = z_f(S.time, pos_1[2])

        delta_pos_0 = np.linalg.norm(S.part[0].pos - pos_0)
        delta_pos_1 = np.linalg.norm(S.part[1].pos - pos_1)

        self.assertLess(1.4e-3, delta_pos_0)
        self.assertLess(delta_pos_0, 1.6e-3)
        self.assertLess(4.9e-4, delta_pos_1)
        self.assertLess(delta_pos_1, 5.1e-4)


if __name__ == '__main__':
    print("Features: ", espressomd.features())
    ut.main()
