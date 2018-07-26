import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize
import unittest as ut

import espressomd.lb
import espressomd.lbboundaries
import espressomd.shapes

AGRID = .25
EXT_FORCE = .2
VISC = .7
DENS = 1.7
TIME_STEP = 0.1


def poiseuille_flow(z, H, ext_force, dyn_visc):
    return ext_force * 1. / (2 * dyn_visc) * (H**2.0 / 4.0 - z**2.0)


class LBPoiseuilleCommon(object):
    lbf = None
    system = espressomd.System(box_l=[10.0, 3.0, 3.0])
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.4 * AGRID

    def prepare(self):
        self.system.actors.clear()
        self.system.actors.add(self.lbf)
        wall_shape1 = espressomd.shapes.Wall(normal=[1, 0, 0], dist=AGRID)
        wall_shape2 = espressomd.shapes.Wall(
            normal=[-1, 0, 0], dist=-(self.system.box_l[0] - AGRID))
        wall1 = espressomd.lbboundaries.LBBoundary(shape=wall_shape1)
        wall2 = espressomd.lbboundaries.LBBoundary(shape=wall_shape2)

        self.system.lbboundaries.add(wall1)
        self.system.lbboundaries.add(wall2)

        mid_indices = [int((self.system.box_l[0] / AGRID) / 2), int(
            (self.system.box_l[1] / AGRID) / 2), int((self.system.box_l[2] / AGRID) / 2)]
        diff = float("inf")
        old_val = self.lbf[mid_indices].velocity[2]
        while diff > 0.01:
            self.system.integrator.run(100)
            new_val = self.lbf[mid_indices].velocity[2]
            diff = abs(new_val - old_val)
            old_val = new_val

    def test_profile(self):
        self.prepare()
        velocities = np.zeros((int(self.system.box_l[0] / AGRID), 2))

        for x in range(velocities.shape[0]):
            v_tmp = []
            for y in range(int(self.system.box_l[1] + 1)):
                for z in range(int(self.system.box_l[2] + 1)):
                    v_tmp.append(self.lbf[x, y, z].velocity[2])
            velocities[x, 1] = np.mean(np.array(v_tmp))
            velocities[x, 0] = x * AGRID

        v_measured = velocities[1:-1, 1]
        v_expected = poiseuille_flow(velocities[1:-1, 0] + 0.5 * AGRID - 0.5 *
                                     self.system.box_l[0], self.system.box_l[0] - 2.0 * AGRID, EXT_FORCE, VISC * DENS)
        rmsd = np.sqrt(np.sum(np.square(v_expected - v_measured)))
        self.assertLess(rmsd, 0.1*AGRID/TIME_STEP)

@ut.skipIf(not espressomd.has_features(['LB', 'LB_BOUNDARIES']), "Skipping test due to missing features.")
class LBCPUPoiseuille(ut.TestCase, LBPoiseuilleCommon):
    lbf = espressomd.lb.LBFluid(agrid=AGRID, dens=DENS, visc=VISC,
                                fric=1.0, tau=TIME_STEP, ext_force_density=[0.0, 0.0, EXT_FORCE])


@ut.skipIf(not espressomd.has_features(['LB_GPU', 'LB_BOUNDARIES_GPU']), "Skipping test due to missing features.")
class LBGPUPoiseuille(ut.TestCase, LBPoiseuilleCommon):
    lbf = espressomd.lb.LBFluidGPU(agrid=AGRID, dens=DENS, visc=VISC,
                                   fric=1.0, tau=TIME_STEP, ext_force_density=[0.0, 0.0, EXT_FORCE])


if __name__ == '__main__':
    ut.main()
