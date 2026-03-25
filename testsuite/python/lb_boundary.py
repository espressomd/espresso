#
# Copyright (C) 2010-2026 The ESPResSo project
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
import espressomd
import espressomd.lb
import espressomd.shapes
import numpy as np


class LBBoundariesBase:
    system = espressomd.System(box_l=[10.0, 5.0, 5.0])
    system.time_step = 0.01
    system.cell_system.skin = 0.1

    wall_shape1 = espressomd.shapes.Wall(normal=[1., 0., 0.], dist=2.5)
    wall_shape2 = espressomd.shapes.Wall(normal=[-1., 0., 0.], dist=-7.5)

    def setUp(self):
        self.lbf = self.lb_class(
            kinematic_viscosity=1.0, density=1.0, agrid=0.5, tau=1.0,
            **self.lb_params)
        self.system.lb = self.lbf
        self.system.thermostat.set_lb(LB_fluid=self.lbf, seed=3, gamma=1.)

    def tearDown(self):
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.lb = None

    def check_boundary_flags(self, slip_velocity1, slip_velocity2):
        def vbb2vel(values):
            velocities = np.empty((*values.shape, 3), dtype=float)
            for index in np.ndindex(*values.shape):
                velocities[index] = values[index].velocity
            return velocities
        lbb1 = self.lbf[:5, :, :]
        lbb2 = self.lbf[15:, :, :]
        lbb3 = self.lbf[5:15, :, :]
        ref_velocity1 = np.tile(slip_velocity1, [5, 10, 10, 1])
        ref_velocity2 = np.tile(slip_velocity2, [5, 10, 10, 1])
        np.testing.assert_equal(np.copy(lbb1.is_boundary), True)
        np.testing.assert_equal(np.copy(lbb2.is_boundary), True)
        np.testing.assert_equal(np.copy(lbb3.is_boundary), False)
        np.testing.assert_allclose(np.copy(lbb1.velocity), ref_velocity1)
        np.testing.assert_allclose(np.copy(lbb2.velocity), ref_velocity2)
        np.testing.assert_allclose(vbb2vel(lbb1.boundary), ref_velocity1)
        np.testing.assert_allclose(vbb2vel(lbb2.boundary), ref_velocity2)
        self.assertTrue(self.lbf[4, 0, 0].is_boundary)
        self.assertFalse(self.lbf[5, 0, 0].is_boundary)
        self.assertFalse(self.lbf[14, 0, 0].is_boundary)
        self.assertTrue(self.lbf[15, 0, 0].is_boundary)
        self.lbf.clear_boundaries()
        np.testing.assert_equal(np.copy(self.lbf[:, :, :].is_boundary), False)

    def test_boundary_flags(self):
        slip_velocity1 = 1e-3 * np.array([1., 2., 3.])
        slip_velocity2 = 1e-3 * np.array([4., 5., 6.])
        value_shape = tuple(self.lbf.shape) + (3,)
        slip_velocity2_all = slip_velocity2 * np.ones(value_shape)
        self.lbf.add_boundary_from_shape(self.wall_shape1, slip_velocity1)
        self.lbf.add_boundary_from_shape(self.wall_shape2, slip_velocity2_all)
        self.check_boundary_flags(slip_velocity1, slip_velocity2)

    def test_union(self):
        union = espressomd.shapes.Union()
        union.add([self.wall_shape1, self.wall_shape2])

        slip_velocity = 1e-3 * np.array([1., 2., 3.])
        self.lbf.add_boundary_from_shape(union, slip_velocity)
        self.check_boundary_flags(slip_velocity, slip_velocity)

    def test_exceptions(self):
        with self.assertRaisesRegex(TypeError, "Parameter 'boundary_type' must be a subclass of VelocityBounceBack"):
            self.lbf.add_boundary_from_shape(
                shape=self.wall_shape1, velocity=[0., 0., 0.],
                boundary_type=self.lb_class)
        with self.assertRaisesRegex(ValueError, "expected an espressomd.shapes.Shape"):
            self.lbf.add_boundary_from_shape(
                shape=self.lbf, velocity=[0., 0., 0.],
                boundary_type=espressomd.lb.VelocityBounceBack)
        with self.assertRaisesRegex(ValueError, r"Cannot process velocity value grid of shape \(4,\)"):
            self.lbf.add_boundary_from_shape(
                shape=self.wall_shape1, velocity=[0., 0., 0., 0.],
                boundary_type=espressomd.lb.VelocityBounceBack)
        self.lbf.add_boundary_from_shape(self.wall_shape1, [0., 0., 0.])

    def test_velocity_interpolation(self):
        xdata = np.linspace(-1, 2, 151)
        slip_velocity = 2e-2 * np.array([1., 2., 3.])

        def get_analytic(axis):
            return np.maximum(0, 1 - np.abs(xdata - 0.5)) * slip_velocity[axis]

        # set up a wall boundary and particles along its normal
        wall = espressomd.shapes.Wall(normal=[1., 0., 0.], dist=self.lbf.agrid)
        self.lbf.add_boundary_from_shape(wall, slip_velocity)
        positions = [[f * self.lbf.agrid, 1., 0.] for f in xdata]
        particles = self.system.part.add(pos=positions)
        if espressomd.has_features("VIRTUAL_SITES_INERTIALESS_TRACERS"):
            tracers = self.system.part.add(pos=positions, propagation=len(
                positions) * [espressomd.propagation.Propagation.TRANS_LB_TRACER])

        # calculate interpolated velocities
        self.system.integrator.run(1)
        lb_vel = []
        part_f = []
        tracer_v = []
        for p in particles:
            lb_vel.append(self.lbf.get_interpolated_velocity(pos=p.pos))
            part_f.append(p.f)
        if espressomd.has_features("VIRTUAL_SITES_INERTIALESS_TRACERS"):
            for p in tracers:
                tracer_v.append(p.v)
        lb_vel = np.copy(lb_vel)
        part_f = np.copy(part_f)
        tracer_v = np.copy(tracer_v)

        rtol = 1e-4 if self.lb_params["single_precision"] else 1e-10
        for i in range(3):
            ref = get_analytic(i)
            np.testing.assert_allclose(lb_vel[:, i], ref, rtol=rtol)
            np.testing.assert_allclose(part_f[:, i], ref, rtol=rtol)
            if espressomd.has_features("VIRTUAL_SITES_INERTIALESS_TRACERS"):
                np.testing.assert_allclose(tracer_v[:, i], ref, rtol=rtol)


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBBoundariesWalberlaDoublePrecisionCPU(LBBoundariesBase, ut.TestCase):
    lb_class = espressomd.lb.LBFluid
    lb_params = {"single_precision": False, "gpu": False}


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBBoundariesWalberlaSinglePrecisionCPU(LBBoundariesBase, ut.TestCase):
    lb_class = espressomd.lb.LBFluid
    lb_params = {"single_precision": True, "gpu": False}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class LBBoundariesWalberlaDoublePrecisionGPU(LBBoundariesBase, ut.TestCase):
    lb_class = espressomd.lb.LBFluid
    lb_params = {"single_precision": False, "gpu": True}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class LBBoundariesWalberlaSinglePrecisionGPU(LBBoundariesBase, ut.TestCase):
    lb_class = espressomd.lb.LBFluid
    lb_params = {"single_precision": True, "gpu": True}


@utx.skipIfMissingFeatures(["WALBERLA"])
class LBBoundariesWalberlaDoublePrecisionCPU(LBBoundariesBase, ut.TestCase):
    lb_class = espressomd.lb.LBFluid
    lb_params = {"single_precision": False, "gpu": False,
                 "blocks_per_mpi_rank": [2, 1, 1]}


if __name__ == "__main__":
    ut.main()
