#
# Copyright (C) 2010-2022 The ESPResSo project
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

import unittest as ut
import unittest_decorators as utx
import numpy as np

import espressomd.lb
import espressomd.shapes

AGRID = .25
EXT_FORCE = .1
KINEMATIC_VISC = 2.7
DENS = 1.7
TIME_STEP = 0.07
LB_PARAMS = {'agrid': AGRID,
             'density': DENS,
             'kinematic_viscosity': KINEMATIC_VISC,
             'tau': TIME_STEP,
             'ext_force_density': [0.0, 0.0, EXT_FORCE]}


def poiseuille_flow(z, H, ext_force_density, dyn_visc):
    """
    Analytical solution for planar Poiseuille flow.

    Parameters
    ----------
    z : :obj:`float`
        Distance to the mid plane of the channel.
    H : :obj:`float`
        Distance between the boundaries.
    ext_force_density : :obj:`float`
        Force density on the fluid normal to the boundaries.
    dyn_visc : :obj:`float`
        Dynamic viscosity of the LB fluid.

    """
    return ext_force_density * 1. / (2 * dyn_visc) * (H**2.0 / 4.0 - z**2.0)


class LBPoiseuilleCommon:

    """
    Check the lattice-Boltzmann pressure-driven flow in a slab system
    by comparing to the analytical solution for the planar Poiseuille.
    """

    system = espressomd.System(box_l=[9.0, 3.0, 3.0])
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.4 * AGRID

    def setUp(self):
        self.lbf = self.lb_class(**LB_PARAMS, **self.lb_params)
        self.system.actors.add(self.lbf)

    def tearDown(self):
        self.system.actors.clear()

    def prepare(self):
        """
        Integrate the LB fluid until steady state is reached within a certain
        accuracy.

        """
        wall_shape1 = espressomd.shapes.Wall(normal=[1, 0, 0], dist=AGRID)
        wall_shape2 = espressomd.shapes.Wall(
            normal=[-1, 0, 0], dist=-(self.system.box_l[0] - AGRID))

        self.lbf.add_boundary_from_shape(wall_shape1)
        self.lbf.add_boundary_from_shape(wall_shape2)

        mid_indices = (self.system.box_l / AGRID / 2).astype(int)
        diff = float("inf")
        old_val = self.lbf[mid_indices].velocity[2]
        while diff > 0.005:
            self.system.integrator.run(200)
            new_val = self.lbf[mid_indices].velocity[2]
            diff = abs(new_val - old_val)
            old_val = new_val

    def test_profile(self):
        """
        Compare against analytical function by calculating the RMSD.

        """
        self.prepare()
        velocities = np.zeros((int(self.system.box_l[0] / AGRID), 2))

        for x in range(velocities.shape[0]):
            v_tmp = []
            for y in range(int(self.system.box_l[1] + 1)):
                for z in range(int(self.system.box_l[2] + 1)):
                    v_tmp.append(self.lbf[x, y, z].velocity[2])
            velocities[x, 1] = np.mean(np.array(v_tmp))
            velocities[x, 0] = (x + 0.5) * AGRID

        v_measured = velocities[1:-1, 1]
        v_expected = poiseuille_flow(velocities[1:-1, 0] - 0.5 * self.system.box_l[0],
                                     self.system.box_l[0] - 2.0 * AGRID,
                                     EXT_FORCE,
                                     KINEMATIC_VISC * DENS)
        np.testing.assert_allclose(v_measured, v_expected, rtol=5E-5)


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBPoiseuilleWalberla(LBPoiseuilleCommon, ut.TestCase):

    """Test for the Walberla implementation of the LB in double-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": False}


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBPoiseuilleWalberlaSinglePrecision(LBPoiseuilleCommon, ut.TestCase):

    """Test for the Walberla implementation of the LB in single-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": True}


if __name__ == '__main__':
    ut.main()
