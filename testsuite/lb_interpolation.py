import unittest as ut
import numpy as np

import espressomd
import espressomd.shapes
import espressomd.lb

AGRID = 1.0
VISC = 1.0
DENS = 1.0
FRIC = 1.0
TAU = 0.1
EXT_FORCE_DENSITY = [0.0, 0.0, 0.1]
BOX_L = 10.0
TIME_STEP = TAU
LB_PARAMETERS = {
    'agrid': AGRID,
    'visc': VISC,
    'dens': DENS,
    'fric': FRIC,
    'tau': TAU,
    'ext_force_density': EXT_FORCE_DENSITY}


class LBInterpolation(object):
    lbf = None
    system = espressomd.System(box_l=[10.0] * 3)
    system.cell_system.skin = 0.4 * AGRID
    system.time_step = TIME_STEP

    def set_boundaries(self):
        """Place boundaries *not* exactly on a LB node.

        """
        wall_shape1 = espressomd.shapes.Wall(
            normal=[1, 0, 0], dist=0.6 * AGRID)
        wall_shape2 = espressomd.shapes.Wall(
            normal=[-1, 0, 0], dist=-(BOX_L - 0.6 * AGRID))
        self.system.lbboundaries.add(
            espressomd.lbboundaries.LBBoundary(shape=wall_shape1))
        self.system.lbboundaries.add(
            espressomd.lbboundaries.LBBoundary(shape=wall_shape2))

    def test_interpolated_velocity(self):
        """
        Check that the interpolated LB fluid velocity is zero between boundary
        node and first fluid node.

        """
        pass


@ut.skipIf(not espressomd.has_features(['LB', 'LB_BOUNDARIES']), "Skipped, features missing.")
class LBInterpolationCPU(ut.TestCase, LBInterpolation):
    def setUp(self):
        self.lbf = espressomd.lb.LBFluid(**LB_PARAMETERS)


@ut.skipIf(not espressomd.has_features(['LB_GPU', 'LB_BOUNDARIES_GPU']), "Skipped, features missing.")
class LBInterpolationGPU(ut.TestCase, LBInterpolation):
    def setUp(self):
        self.lbf = espressomd.lb.LBFluidGPU(**LB_PARAMETERS)


if __name__ == "__main__":
    ut.main()
