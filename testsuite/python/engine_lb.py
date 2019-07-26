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
import unittest_decorators as utx
from tests_common import abspath

import espressomd
import espressomd.lb

import os
import numpy as np

try:
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy
except ImportError:
    pass


def calculate_vtk_max_pointwise_difference(file1, file2, tol=1e-6):
    arrays = [0] * 2

    reader = vtk.vtkStructuredPointsReader()
    for i, fname in enumerate([file1, file2]):
        reader.SetFileName(fname)
        reader.Update()
        data = reader.GetOutput().GetPointData()
        arrays[i] = np.array([vtk_to_numpy(data.GetArray(n))
                              for n in range(data.GetNumberOfArrays())])

    try:
        return np.allclose(arrays[0], arrays[1], rtol=0, atol=tol), np.max(
            np.abs(arrays[0] - arrays[1]))
    except BaseException:
        return False, np.inf


@utx.skipIfMissingFeatures(["ENGINE"])
@utx.skipIfMissingModules(['vtk'])
class SwimmerTest(ut.TestCase):
    S = espressomd.System(box_l=[1.0, 1.0, 1.0])
    S.seed = S.cell_system.get_state()['n_nodes'] * [1234]

    def setUp(self):
        S = self.S
        S.box_l = 3 * [12]
        S.cell_system.skin = 0.1
        S.time_step = 0.01

        S.part.add(id=0, pos=[6.0, 3.0, 2.0], quat=np.sqrt([.5, .5, 0, 0]),
                   swimming={"mode": "pusher", "v_swim": 0.10,
                             "dipole_length": 1.0, "rotational_friction": 2.0})
        S.part.add(id=1, pos=[2.0, 3.0, 6.0], quat=np.sqrt([.5, 0, .5, 0]),
                   swimming={"mode": "pusher", "f_swim": 0.03,
                             "dipole_length": 2.0, "rotational_friction": 20.})
        S.part.add(id=2, pos=[3.0, 2.0, 6.0], quat=np.sqrt([.5, 0, 0, .5]),
                   swimming={"mode": "puller", "v_swim": 0.15,
                             "dipole_length": 0.5, "rotational_friction": 15.})
        S.part.add(id=3, pos=[3.0, 6.0, 2.0], quat=np.sqrt([0, 0, .5, .5]),
                   swimming={"mode": "puller", "f_swim": 0.05,
                             "dipole_length": 1.5, "rotational_friction": 6.0})
        S.part[:].rotation = [1, 1, 1]

    def tearDown(self):
        self.S.part.clear()

    def run_and_check(self, lbm, vtk_ref, vtk_out, tol):
        self.S.integrator.run(2000)

        lbm.print_vtk_velocity(vtk_out)
        identical, err_max = calculate_vtk_max_pointwise_difference(
            vtk_ref, vtk_out, tol=tol)
        os.remove(vtk_out)
        print(
            "Maximum deviation to the reference point is: {}".format(err_max))
        self.assertTrue(identical)


class SwimmerTestCPU(SwimmerTest):

    def test(self):
        lbm = espressomd.lb.LBFluid(
            agrid=1.0, tau=self.S.time_step, visc=1.0, dens=1.0)
        self.S.actors.add(lbm)
        self.S.thermostat.set_lb(LB_fluid=lbm, gamma=0.5)
        self.run_and_check(
            lbm, abspath("data/engine_lb.vtk"),
            "engine_test_cpu_tmp.vtk", 1.5e-6)
        self.S.thermostat.turn_off()
        self.S.actors.remove(lbm)


@utx.skipIfMissingGPU()
class SwimmerTestGPU(SwimmerTest):

    def test(self):
        lbm = espressomd.lb.LBFluidGPU(
            agrid=1.0,
            tau=self.S.time_step,
            visc=1.0,
            dens=1.0
        )
        self.S.actors.add(lbm)
        self.S.thermostat.set_lb(LB_fluid=lbm, gamma=0.5)
        self.run_and_check(
            lbm, abspath("data/engine_lbgpu_2pt.vtk"),
            "engine_test_gpu_tmp.vtk", 2.0e-7)
        self.S.thermostat.turn_off()
        self.S.actors.remove(lbm)


if __name__ == '__main__':
    ut.main()
