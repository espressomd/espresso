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

AGRID = 0.6
VISC = 5.2
DENS = 2.3
TIME_STEP = 0.02
# Box size will be H +2 AGRID to make room for walls.
# The number of grid cells should be divisible by four and 3 in all directions
# for testing on multiple mpi nodes.
H = 12 * AGRID
W = 6 * AGRID
SHEAR_VELOCITY = 0.3

LB_PARAMS = {'agrid': AGRID,
             'density': DENS,
             'viscosity': VISC,
             'tau': TIME_STEP
             }


def shear_flow(x, t, nu, v, h, k_max):
    """
    Analytical solution for driven shear flow between two plates.

    Parameters
    ----------
    x : :obj:`float`
        Position from the left plane.
    t : :obj:`float`
        Time since start of the shearing.
    nu : :obj:`float`
        Kinematic viscosity.
    v : :obj:`float`
        Shear rate.
    h : :obj:`float`
        Distance between the plates.
    k_max : :obj:`int`
        Maximum considered wave number.

    Returns
    -------
    :obj:`double` : Analytical velocity

    """

    u = x / h - 0.5
    for k in np.arange(1, k_max + 1):
        u += 1.0 / (np.pi * k) * np.exp(
            -4 * np.pi ** 2 * nu * k ** 2 / h ** 2 * t) * np.sin(2 * np.pi / h * k * x)
    return -v * u


class LBShearCommon:

    """
    Check the lattice-Boltzmann lid-driven shear flow in a slab system
    by comparing to the analytical solution.
    """
    system = espressomd.System(box_l=[H + 2. * AGRID, W, W])
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.4 * AGRID

    def setUp(self):
        self.lbf = self.lb_class(**LB_PARAMS, **self.lb_params)

    def tearDown(self):
        self.system.actors.clear()

    def check_profile(self, shear_plane_normal, shear_direction):
        """
        Integrate the LB fluid and regularly compare with
        the exact solution.

        """
        self.tearDown()
        self.setUp()
        self.system.box_l = np.max(
            ((W, W, W), shear_plane_normal * (H + 2 * AGRID)), 0)
        self.system.actors.add(self.lbf)
        self.lbf.clear_boundaries()

        wall_shape1 = espressomd.shapes.Wall(
            normal=shear_plane_normal, dist=AGRID)
        wall_shape2 = espressomd.shapes.Wall(
            normal=-1.0 * shear_plane_normal, dist=-(H + AGRID))

        self.lbf.add_boundary_from_shape(
            wall_shape1, velocity=-.5 * SHEAR_VELOCITY * shear_direction)
        self.lbf.add_boundary_from_shape(
            wall_shape2, velocity=.5 * SHEAR_VELOCITY * shear_direction)

        t0 = self.system.time
        sample_points = int(H / AGRID - 1)

        # warmup
        self.system.integrator.run(40)
        for _ in range(9):
            self.system.integrator.run(20)

            shear_profile = shear_flow(
                x=(np.arange(0, sample_points) + .5) * AGRID,
                t=self.system.time - t0,
                nu=VISC,
                v=SHEAR_VELOCITY,
                h=H,
                k_max=100)
            v_expected = -np.outer(shear_profile, shear_direction)
            v_measured = np.zeros((sample_points, 3))
            for j in range(0, sample_points):
                ind = np.max(((1, 1, 1), shear_plane_normal * j + 1), 0)
                v_measured[j] = self.lbf[ind[0], ind[1], ind[2]].velocity
            np.testing.assert_allclose(v_measured, v_expected, atol=8E-4)

        # speed of sound of the LB fluid in MD units (agrid/tau is due to
        # LB->MD unit conversion)
        speed_of_sound = 1. / np.sqrt(3.) * self.lbf.agrid / self.lbf.tau
        # equation of state for the LB fluid
        p_eq = speed_of_sound**2 * DENS
        # see Eq. 1.15 and 1.29 in Timm Krüger et al. "The lattice Boltzmann
        # method", Springer International Publishing 10 (2017): 978-3,
        # https://doi.org/10.1007/978-3-319-44649-3
        # and
        # https://de.wikipedia.org/wiki/Navier-Stokes-Gleichungen
        # note that for an incompressible fluid the viscous stress tensor is
        # defined as \sigma = -p 1 + \mu [\nabla * u + (\nabla * u)^T]
        # where 'p' is the static pressure, '\mu' is the dynamic viscosity,
        # '*' denotes the outer product and 'u' is the velocity field
        # NOTE: the so called stress property of the fluid is actually the
        # pressure tensor not the viscous stress tensor!
        shear_rate = SHEAR_VELOCITY / H
        dynamic_viscosity = self.lbf.viscosity * DENS
        p_expected = p_eq * np.identity(3) - dynamic_viscosity * shear_rate * (
            np.outer(shear_plane_normal, shear_direction) + np.transpose(np.outer(shear_plane_normal, shear_direction)))
        # Walberla TODO
        for n in []:  # (2, 3, 4), (3, 4, 2), (5, 4, 3):
            node_pressure_tensor = np.copy(
                self.lbf[n[0], n[1], n[2]].pressure_tensor)
            np.testing.assert_allclose(node_pressure_tensor,
                                       p_expected, atol=4E-5, rtol=5E-3)  # WALBERLA TODO

        # TODO: WALBERLA: (#4381) boundary forces not reliable at the moment
#        np.testing.assert_allclose(
#            np.copy(wall1.get_force()),
#            -np.copy(wall2.get_force()),
#            atol=1E-4)
#        np.testing.assert_allclose(np.dot(np.copy(wall1.get_force()), shear_direction),
# SHEAR_VELOCITY / H * W**2 * dynamic_viscosity, atol=2E-4)

    def test(self):
        x = np.array((1, 0, 0), dtype=int)
        y = np.array((0, 1, 0), dtype=int)
        z = np.array((0, 0, 1), dtype=int)
        self.check_profile(x, y)
        self.check_profile(x, z)
        self.check_profile(y, z)
        self.check_profile(x, -y)
        self.check_profile(x, -z)
        self.check_profile(y, -z)


@utx.skipIfMissingFeatures(['WALBERLA'])
class LBShearWalberla(LBShearCommon, ut.TestCase):

    """Test for the Walberla implementation of the LB in double-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': False}


@utx.skipIfMissingFeatures(['WALBERLA'])
class LBShearWalberlaSinglePrecision(LBShearCommon, ut.TestCase):

    """Test for the Walberla implementation of the LB in single-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': True}


if __name__ == '__main__':
    ut.main()
