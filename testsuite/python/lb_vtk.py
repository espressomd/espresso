#
# Copyright (C) 2010-2020 The ESPResSo project
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

import os
import numpy as np

try:
    import vtk
    from vtk.util import numpy_support as VN
    skipIfMissingPythonPackage = ut.case._id
except ImportError:
    skipIfMissingPythonPackage = ut.skip(
        "Python module vtk not available, skipping test!")

import espressomd
import espressomd.lb
if espressomd.has_features('LB_BOUNDARIES'):
    import espressomd.lbboundaries
    import espressomd.shapes


@skipIfMissingPythonPackage
class TestVTK:
    system = espressomd.System(box_l=3 * [16])
    system.time_step = 0.01
    system.cell_system.skin = 0.4

    def tearDown(self):
        self.system.actors.clear()
        self.system.thermostat.turn_off()

    def parse_vtk(self, filepath, name, shape):
        reader = vtk.vtkStructuredPointsReader()
        reader.SetFileName(filepath)
        reader.ReadAllVectorsOn()
        reader.ReadAllScalarsOn()
        reader.Update()

        data = reader.GetOutput()
        points = data.GetPointData()

        return VN.vtk_to_numpy(points.GetArray(name)).reshape(shape, order='F')

    def test_vtk(self):
        '''
        Check VTK files.
        '''

        # setup LB system
        lbf = self.lb_class(
            kT=0, agrid=1.0, dens=1.0, visc=1.0, tau=0.1,
            ext_force_density=[0, 0.03, 0])
        self.system.actors.add(lbf)
        shape = [16, 16, 16]
        if espressomd.has_features('LB_BOUNDARIES'):
            self.system.lbboundaries.add(espressomd.lbboundaries.LBBoundary(
                shape=espressomd.shapes.Wall(normal=[1, 0, 0], dist=1.5)))
            self.system.lbboundaries.add(espressomd.lbboundaries.LBBoundary(
                shape=espressomd.shapes.Wall(normal=[-1, 0, 0], dist=-14.5)))

        os.makedirs('vtk_out', exist_ok=True)
        filepaths = ['vtk_out/boundary.vtk', 'vtk_out/velocity.vtk',
                     'vtk_out/velocity_bb.vtk', 'vtk_out/delme']

        # cleanup action
        for filepath in filepaths:
            if os.path.exists(filepath):
                os.remove(filepath)

        # write VTK files
        self.system.integrator.run(100)
        with self.assertRaises(RuntimeError):
            lbf.print_vtk_velocity('non_existent_folder/file')
        lbf.print_vtk_boundary('vtk_out/boundary.vtk')
        lbf.print_vtk_velocity('vtk_out/velocity.vtk')
        with self.assertRaises(ValueError):
            lbf.print_vtk_velocity('vtk_out/delme', 3 * [0], None)
        with self.assertRaises(ValueError):
            lbf.print_vtk_velocity('vtk_out/delme', None, 3 * [0])
        with self.assertRaises(RuntimeError):
            lbf.print_vtk_velocity('vtk_out/delme', [-2, 1, 1], 3 * [1])
        with self.assertRaises(RuntimeError):
            lbf.print_vtk_velocity('vtk_out/delme', 3 * [0], [1, 2, 16])
        with self.assertRaises(ValueError):
            lbf.print_vtk_velocity('vtk_out/delme', [1, 1], 3 * [1])
        with self.assertRaises(ValueError):
            lbf.print_vtk_velocity('vtk_out/delme', 3 * [1], np.array([2, 3]))
        bb1, bb2 = ([1, 2, 3], [12, 13, 14])
        lbf.print_vtk_velocity('vtk_out/velocity_bb.vtk', bb1, bb2)

        # check VTK files exist
        for filepath in filepaths:
            self.assertTrue(
                os.path.exists(filepath),
                'VTK file "{}" not written to disk'.format(filepath))

        # check VTK values match node values
        node_velocity = np.zeros(shape + [3])
        node_boundary = np.zeros(shape, dtype=int)
        for i in range(shape[0]):
            for j in range(shape[1]):
                for k in range(shape[2]):
                    node = lbf[i, j, k]
                    node_velocity[i, j, k] = node.velocity
                    node_boundary[i, j, k] = node.boundary
        node_velocity_bb = node_velocity[bb1[0]:bb2[0] + 1,
                                         bb1[1]:bb2[1] + 1,
                                         bb1[2]:bb2[2] + 1]

        vtk_velocity = self.parse_vtk('vtk_out/velocity.vtk', 'velocity',
                                      node_velocity.shape)
        np.testing.assert_allclose(vtk_velocity, node_velocity, atol=5e-7)

        vtk_velocity_bb = self.parse_vtk('vtk_out/velocity_bb.vtk', 'velocity',
                                         node_velocity_bb.shape)
        np.testing.assert_allclose(
            vtk_velocity_bb, node_velocity_bb, atol=5e-7)

        vtk_boundary = self.parse_vtk(
            'vtk_out/boundary.vtk', 'boundary', shape)
        np.testing.assert_equal(vtk_boundary, node_boundary.astype(int))


class TestVTKCPU(TestVTK, ut.TestCase):

    def setUp(self):
        self.lb_class = espressomd.lb.LBFluid


@utx.skipIfMissingGPU()
class TestVTKGPU(TestVTK, ut.TestCase):

    def setUp(self):
        self.lb_class = espressomd.lb.LBFluidGPU


if __name__ == '__main__':
    ut.main()
