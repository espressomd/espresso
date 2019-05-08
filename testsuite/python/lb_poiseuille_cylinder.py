# Copyright (C) 2010-2018 The ESPResSo project
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
import unittest as ut
import numpy as np

import espressomd.lb
import espressomd.lbboundaries
import espressomd.shapes


"""
Check the Lattice Boltzmann 'pressure' driven flow in a cylindrical constraint
by comparing to the analytical solution.


"""


AGRID = .25
EXT_FORCE = .1
VISC = 2.7
DENS = 1.7
TIME_STEP = 0.1
LB_PARAMS = {'agrid': AGRID,
             'dens': DENS,
             'visc': VISC,
             'tau': TIME_STEP,
             'ext_force_density': [0.0, 0.0, EXT_FORCE]}


def poiseuille_flow(r, R, ext_force_density, dyn_visc):
    """
    Analytical solution for Poiseuille flow.

    Parameters
    ----------
    r : :obj:`float`
        Distance to the center of the tube.
    R : :obj:`float`
        Radius of the tube.
    ext_force_density : :obj:`float`
        Force density on the fluid parallel to the boundaries.
    dyn_visc : :obj:`float`
        Dynamic viscosity of the fluid.

    """
    return ext_force_density * 1. / (4 * dyn_visc) * (R**2.0 - r**2.0)


class LBPoiseuilleCommon(object):

    """Base class of the test that holds the test logic."""
    lbf = None
    system = espressomd.System(box_l=[9.0, 9.0, 9.0])
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.4 * AGRID

    def prepare(self):
        """
        Integrate the LB fluid until steady state is reached within a certain
        accuracy.

        """
        self.system.actors.clear()
        self.system.actors.add(self.lbf)
        cylinder_shape = espressomd.shapes.Cylinder(center=self.system.box_l / 2.0, axis=[0, 0, 1], direction=-1, radius=self.system.box_l[2] / 2.0 - 1.0, length=self.system.box_l[2] * 1.5)
        cylinder = espressomd.lbboundaries.LBBoundary(shape=cylinder_shape)

        self.system.lbboundaries.add(cylinder)

        mid_indices = [int((self.system.box_l[0] / AGRID) / 2),
                       int((self.system.box_l[1] / AGRID) / 2),
                       int((self.system.box_l[2] / AGRID) / 2)]
        diff = float("inf")
        old_val = self.lbf[mid_indices].velocity[2]
        while diff > 0.001:
            self.system.integrator.run(200)
            new_val = self.lbf[mid_indices].velocity[2]
            diff = abs(new_val - old_val)
            old_val = new_val

    def test_profile(self):
        """
        Compare against analytical function by calculating the RMSD.

        """
        self.prepare()
        velocities = np.zeros(int(self.system.box_l[0] / AGRID))
        positions = np.zeros_like(velocities)

        for y in range(velocities.shape[0]):
            v_tmp = []
            for z in range(int(self.system.box_l[2] / AGRID)):
                v_tmp.append(self.lbf[int(self.system.box_l[0] / AGRID) / 2, y, z].velocity[2])
            velocities[y] = np.mean(np.array(v_tmp))
            positions[y] = (y + 0.5) * AGRID

        v_measured = velocities[1:-1]
        v_expected = poiseuille_flow(positions[1:-1] - 0.5 * self.system.box_l[0],
                                     self.system.box_l[0] / 2.0 - 1.0,
                                     EXT_FORCE,
                                     VISC * DENS)
        rmsd = np.sqrt(np.sum(np.square(v_expected - v_measured)))
        self.assertLess(rmsd, 0.02 * AGRID / TIME_STEP)


@ut.skipIf(not espressomd.has_features(
    ['LB_BOUNDARIES', 'EXTERNAL_FORCES']), "Skipping test due to missing features.")
class LBCPUPoiseuille(ut.TestCase, LBPoiseuilleCommon):

    """Test for the CPU implementation of the LB."""

    def setUp(self):
        self.lbf = espressomd.lb.LBFluid(**LB_PARAMS)


@ut.skipIf(not espressomd.gpu_available() or not espressomd.has_features(
    ['CUDA', 'LB_BOUNDARIES_GPU', 'EXTERNAL_FORCES']), "Skipping test due to missing features or gpu.")
class LBGPUPoiseuille(ut.TestCase, LBPoiseuilleCommon):

    """Test for the GPU implementation of the LB."""

    def setUp(self):
        self.lbf = espressomd.lb.LBFluidGPU(**LB_PARAMS)


if __name__ == '__main__':
    ut.main()
