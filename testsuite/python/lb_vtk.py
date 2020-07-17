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
import vtk
from vtk.util import numpy_support as VN

import espressomd
import espressomd.lb
if espressomd.has_features('LB_BOUNDARIES'):
    import espressomd.lbboundaries
    import espressomd.shapes


@utx.skipIfMissingFeatures("LB_WALBERLA")
class TestVTK(ut.TestCase):
    system = espressomd.System(box_l=3 * [16])
    system.time_step = 0.01
    system.cell_system.skin = 0.4
    lbf = espressomd.lb.LBFluidWalberla(kT=0, agrid=1.0, dens=1.0, visc=1.0,
                                        tau=0.1, ext_force_density=[0, 0.03, 0])

    def get_cell_array(self, cell, name, shape):
        return VN.vtk_to_numpy(cell.GetArray(name)).reshape(shape, order='F')

    def test_vtk(self):
        filepath = 'vtk_out/test_lb_vtk/simulation_step_9.vtu'

        # setup LB system
        self.system.actors.add(self.lbf)
        if espressomd.has_features('LB_BOUNDARIES'):
            self.system.lbboundaries.add(espressomd.lbboundaries.LBBoundary(
                shape=espressomd.shapes.Wall(normal=[1, 0, 0], dist=1.5)))
            self.system.lbboundaries.add(espressomd.lbboundaries.LBBoundary(
                shape=espressomd.shapes.Wall(normal=[-1, 0, 0], dist=14.5)))

        # cleanup action
        if os.path.exists(filepath):
            os.remove(filepath)

        # write VTK file
        vtk_observables = ['density', 'velocity_vector', 'pressure_tensor']
        self.lbf.write_vtk('test_lb_vtk', vtk_observables, 1)
        self.system.integrator.run(100)

        # parse VTK file
        self.assertTrue(
            os.path.exists(filepath),
            'VTK file not written to disk')
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(filepath)
        reader.Update()

        data = reader.GetOutput()
        cell = data.GetCellData()

        x_offset = 0
        shape = [16, 16, 16]
        if espressomd.has_features('LB_BOUNDARIES'):
            x_offset = 2
            shape[0] = 12

        vtk_density = self.get_cell_array(cell, 'DensityFromPDF', shape)
        # TODO Walberla
        # vtk_velocity = self.get_cell_array(
        #     cell, 'VelocityFromPDF', shape + [3])
        # vtk_pressure = self.get_cell_array(
        #     cell, 'PressureTensorFromPDF', shape + [3, 3])

        node_density = np.zeros(shape)
        # TODO Walberla
        # node_velocity = np.zeros(shape + [3])
        # node_pressure = np.zeros(shape + [3, 3])
        for i in range(shape[0]):
            for j in range(shape[1]):
                for k in range(shape[2]):
                    node = self.lbf[i + x_offset, j, k]
                    node_density[i, j, k] = node.density
                    # TODO Walberla
                    # node_velocity[i, j, k] = node.velocity
                    # node_pressure[i, j, k] = node.pressure_tensor

        # VTK files are written with float precision
        np.testing.assert_allclose(node_density, vtk_density, rtol=1e-5)
        # TODO Walberla
        # np.testing.assert_allclose(node_velocity, vtk_velocity, rtol=1e-5)
        # np.testing.assert_allclose(node_pressure, vtk_pressure, rtol=1e-5)


if __name__ == '__main__':
    ut.main()
