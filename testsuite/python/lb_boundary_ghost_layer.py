#
# Copyright (C) 2024 The ESPResSo project
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
import scipy.optimize

import espressomd.lb
import espressomd.shapes

AGRID = 0.5
KINEMATIC_VISC = 2.7
DENS = 1.7
TIME_STEP = 0.07
LB_PARAMS = {"agrid": AGRID, "tau": TIME_STEP, "density": DENS,
             "kinematic_viscosity": KINEMATIC_VISC}


class TestCommon:

    system = espressomd.System(box_l=[16.0, 1.0, 1.0])
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.4 * AGRID
    if espressomd.gpu_available():
        system.cuda_init_handle.call_method("set_device_id_per_rank")
    n_nodes = system.cell_system.get_state()["n_nodes"]

    def setUp(self):
        self.lbf = self.lb_class(**LB_PARAMS, **self.lb_params)
        self.system.lb = self.lbf
        self.ubb = espressomd.lb.VelocityBounceBack([0., 0., 1e-5])

    def tearDown(self):
        self.system.lb = None

    def get_profile(self):
        xdata = np.arange(1, self.lbf.shape[0])
        ydata = []
        for x in xdata:
            ydata.append(np.mean(self.lbf[x, :, :].velocity[:, :, 2]))
        return xdata, np.array(ydata)

    def check_profile(self):
        def quadratic(x, a, b, c):
            return a * x**2 + b * x + c

        self.system.integrator.run(40)
        xdata, ydata = self.get_profile()
        popt_ref = (4e-8, -1e-6, 1e-5)
        popt, _ = scipy.optimize.curve_fit(
            quadratic, xdata, ydata, p0=popt_ref)
        rtol = 0.3 if self.lbf.single_precision else 0.1
        np.testing.assert_allclose(popt, popt_ref, rtol=0.5, atol=0.)
        np.testing.assert_allclose(ydata, quadratic(xdata, *popt),
                                   rtol=rtol, atol=0.)

    def test_node_setter(self):
        for i in (0, 1):
            for j in (0, 1):
                self.lbf[0, i, j].boundary = self.ubb
        self.check_profile()

    def test_slice_setter(self):
        self.lbf[0, :, :].boundary = self.ubb
        self.check_profile()

    def test_shape_setter(self):
        shape = espressomd.shapes.Wall(normal=[1, 0, 0], dist=AGRID)
        self.lbf.add_boundary_from_shape(shape, velocity=self.ubb.velocity)
        self.check_profile()


@utx.skipIfMissingFeatures(["WALBERLA"])
@ut.skipIf(TestCommon.n_nodes != 2, "only runs for 2 MPI ranks")
class LBPoiseuilleWalberlaSinglePrecisionCPU(TestCommon, ut.TestCase):
    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": True}


@utx.skipIfMissingFeatures(["WALBERLA"])
@ut.skipIf(TestCommon.n_nodes != 2, "only runs for 2 MPI ranks")
class LBPoiseuilleWalberlaDoublePrecisionCPU(TestCommon, ut.TestCase):
    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": False}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
@ut.skipIf(TestCommon.n_nodes != 2, "only runs for 2 MPI ranks")
class LBPoiseuilleWalberlaSinglePrecisionGPU(TestCommon, ut.TestCase):
    lb_class = espressomd.lb.LBFluidWalberlaGPU
    lb_params = {"single_precision": True}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
@ut.skipIf(TestCommon.n_nodes != 2, "only runs for 2 MPI ranks")
class LBPoiseuilleWalberlaDoublePrecisionGPU(TestCommon, ut.TestCase):
    lb_class = espressomd.lb.LBFluidWalberlaGPU
    lb_params = {"single_precision": False}


if __name__ == "__main__":
    ut.main()
