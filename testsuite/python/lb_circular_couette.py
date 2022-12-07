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
import espressomd.observables
import espressomd.math


AGRID = .5
TIME_STEP = 0.1
GRID_SIZE = np.array([63, 63, 4])


def taylor_couette(v1, v2, r1, r2):
    # Taylor-Couette equation
    omega1 = v1 / r1
    omega2 = v2 / r2
    eta = r1 / r2
    a = (omega2 - omega1 * eta**2) / (1. - eta**2)
    b = r1**2 * (omega1 - omega2) / (1. - eta**2)
    return a, b


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBCircularCouetteCommon:

    system = espressomd.System(box_l=(GRID_SIZE + [1, 1, 0]) * AGRID)
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.1
    system.periodicity = [False, False, True]

    def tearDown(self):
        self.system.actors.clear()

    def test_taylor_couette_flow(self):
        """
        Rotate a rod in a cavity filled with fluid. The solution to the
        Navier-Stokes equation, assuming an infinite rod, is the
        Taylor-Couette equation.
        """

        system = self.system
        lb_fluid = espressomd.lb.LBFluidWalberla(
            agrid=AGRID, density=0.5, viscosity=3.2, tau=system.time_step)
        system.actors.add(lb_fluid)

        # set up two cylinders
        cyl_center = AGRID * (GRID_SIZE // 2 + 0.5) * [1, 1, 0]
        cyl1 = espressomd.shapes.Cylinder(
            center=cyl_center, axis=[0, 0, 1], length=3 * system.box_l[2],
            radius=8.1 * AGRID, direction=1)
        cyl2 = espressomd.shapes.Cylinder(
            center=cyl_center, axis=[0, 0, 1], length=3 * system.box_l[2],
            radius=30.1 * AGRID, direction=-1)
        lb_fluid.add_boundary_from_shape(cyl1)
        lb_fluid.add_boundary_from_shape(cyl2)

        # the system needs to be fully symmetric
        mask = np.copy(lb_fluid[:63, :63, :].is_boundary.astype(int))
        np.testing.assert_array_equal(mask, np.flip(mask, axis=0))
        np.testing.assert_array_equal(mask, np.flip(mask, axis=1))
        np.testing.assert_array_equal(mask, np.flip(mask, axis=2))

        # the system needs to be closed in the x and y directions
        np.testing.assert_array_equal(mask[0, :, :], 1)
        np.testing.assert_array_equal(mask[-1, :, :], 1)
        np.testing.assert_array_equal(mask[:, 0, :], 1)
        np.testing.assert_array_equal(mask[:, -1, :], 1)

        # add tangential slip velocity to the inner cylinder
        slip_vel = 0.01
        surface_nodes = espressomd.lb.edge_detection(
            lb_fluid.get_shape_bitmask(cyl1), system.periodicity)
        tangents = espressomd.lb.calc_cylinder_tangential_vectors(
            cyl1.center, AGRID, 0.5, surface_nodes)
        for node, tangent in zip(surface_nodes, tangents):
            lb_fluid[node].boundary = espressomd.lb.VelocityBounceBack(
                slip_vel * tangent)

        # add observable for the fluid velocity in cylindrical coordinates
        cyl_transform_params = espressomd.math.CylindricalTransformationParameters(
            center=cyl_center, axis=[0, 0, 1], orientation=[1, 0, 0])
        observable = espressomd.observables.CylindricalLBVelocityProfile(
            transform_params=cyl_transform_params,
            n_r_bins=GRID_SIZE[0] // 2,
            n_phi_bins=1,
            n_z_bins=1,
            min_r=0.0,
            max_r=system.box_l[0] / 2,
            min_phi=0.,
            max_phi=2 * np.pi,
            min_z=0.,
            max_z=+system.box_l[2],
            axis=[0.0, 0.0, 1.0],
            sampling_density=1
        )

        # equilibrate the fluid and sample velocities
        obs_data_baseline = observable.calculate()
        system.integrator.run(200)
        obs_data = observable.calculate()
        profile_r = np.copy(observable.bin_centers()).reshape([-1, 3])[:, 0]
        profile_v = np.copy(obs_data - obs_data_baseline).reshape([-1, 3])
        v_r, v_phi, v_z = profile_v.T

        # check velocity is zero for the radial and axial components
        np.testing.assert_allclose(v_r, 0., atol=1e-4)
        np.testing.assert_allclose(v_z, 0., atol=1e-6)

        # check azimuthal velocity is zero inside boundary
        np.testing.assert_allclose(v_phi[:7], 0., atol=1e-7)

        # check azimuthal velocity in the linear regime
        self.assertGreater(v_phi[7], v_phi[6])
        self.assertGreater(v_phi[8], v_phi[7])
        self.assertGreater(v_phi[9], v_phi[8])

        # check azimuthal velocity in the Couette regime
        r = profile_r[10:-1]
        v_phi = v_phi[10:-1]
        a_ref, b_ref = taylor_couette(slip_vel, 0.0, cyl1.radius, cyl2.radius)
        v_phi_ref = a_ref * r + b_ref / r
        v_phi_drift = np.mean(v_phi) - np.mean(v_phi_ref)
        np.testing.assert_allclose(v_phi_drift, 0., atol=1.2e-4)
        np.testing.assert_allclose(v_phi - v_phi_drift, v_phi_ref, atol=1e-4)


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBCircularCouetteWalberla(LBCircularCouetteCommon, ut.TestCase):

    """Test for the Walberla implementation of the LB in double-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': False}


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBCircularCouetteWalberlaSinglePrecision(
        LBCircularCouetteCommon, ut.TestCase):

    """Test for the Walberla implementation of the LB in single-precision."""

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': True}


if __name__ == "__main__":
    ut.main()
