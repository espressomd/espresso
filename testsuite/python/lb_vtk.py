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
    skipIfMissingPythonPackage = utx.no_skip
except ImportError:
    skipIfMissingPythonPackage = ut.skip(
        "Python module vtk not available, skipping test!")

import espressomd
import espressomd.lb
if espressomd.has_features('LB_BOUNDARIES'):
    import espressomd.lbboundaries
    import espressomd.shapes


@skipIfMissingPythonPackage
@utx.skipIfMissingFeatures("LB_WALBERLA")
class TestVTK(ut.TestCase):
    system = espressomd.System(box_l=3 * [16])
    system.time_step = 0.01
    system.cell_system.skin = 0.4

    def get_cell_array(self, cell, name, shape):
        return VN.vtk_to_numpy(cell.GetArray(name)).reshape(shape, order='F')

    def parse_vtk(self, filepath, shape):
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(filepath)
        reader.Update()

        data = reader.GetOutput()
        cell = data.GetCellData()

        vtk_density = self.get_cell_array(cell, 'DensityFromPDF', shape)
        vtk_velocity = self.get_cell_array(
            cell, 'VelocityFromPDF', shape + [3])
        # TODO Walberla
        # vtk_pressure = self.get_cell_array(
        #     cell, 'PressureTensorFromPDF', shape + [3, 3])
        return vtk_density, vtk_velocity

    def test_vtk(self):
        '''
        Check VTK files. Keep in mind VTK files are written with
        float precision.
        '''

        # setup LB system
        self.lbf = espressomd.lb.LBFluidWalberla(
            kT=0, agrid=1.0, dens=1.0, visc=1.0, tau=0.1,
            ext_force_density=[0, 0.03, 0])
        self.system.actors.add(self.lbf)
        x_offset = 0
        shape = [16, 16, 16]
        if espressomd.has_features('LB_BOUNDARIES'):
            self.system.lbboundaries.add(espressomd.lbboundaries.LBBoundary(
                shape=espressomd.shapes.Wall(normal=[1, 0, 0], dist=1.5)))
            self.system.lbboundaries.add(espressomd.lbboundaries.LBBoundary(
                shape=espressomd.shapes.Wall(normal=[-1, 0, 0], dist=-14.5)))
            x_offset = 2
            shape[0] = 12

        n_steps = 100
        lb_steps = int(np.floor(n_steps * self.lbf.tau))
        filepaths = ['vtk_out/test_lb_vtk_end/simulation_step_0.vtu'] + \
            ['vtk_out/test_lb_vtk_continuous/simulation_step_{}.vtu'.format(
                i) for i in range(lb_steps)]

        # cleanup action
        for filepath in filepaths:
            if os.path.exists(filepath):
                os.remove(filepath)

        # write VTK files
        # TODO Walberla: add 'pressure_tensor'
        vtk_obs = ['density', 'velocity_vector']
        self.lbf.add_vtk_writer('test_lb_vtk_continuous', vtk_obs, 1)
        self.system.integrator.run(n_steps)
        lb_vtk = self.lbf.add_vtk_writer('test_lb_vtk_end', vtk_obs, 0)
        lb_vtk.write()

        # check VTK files exist
        for filepath in filepaths:
            self.assertTrue(
                os.path.exists(filepath),
                'VTK file "{}" not written to disk'.format(filepath))

        # check velocity profile is symmetric at all time steps
        for filepath in filepaths:
            vtk_velocity = self.parse_vtk(filepath, shape)[1]
            v_profile = np.mean(
                np.linalg.norm(vtk_velocity, axis=-1),
                axis=(1, 2))
            np.testing.assert_allclose(v_profile, v_profile[::-1], rtol=1e-5)

        # check VTK values match node values in the final time step
        node_density = np.zeros(shape)
        node_velocity = np.zeros(shape + [3])
        # TODO Walberla
        # node_pressure = np.zeros(shape + [3, 3])
        for i in range(shape[0]):
            for j in range(shape[1]):
                for k in range(shape[2]):
                    node = self.lbf[i + x_offset, j, k]
                    node_density[i, j, k] = node.density
                    node_velocity[i, j, k] = node.velocity
                    # TODO Walberla
                    # node_pressure[i, j, k] = node.pressure_tensor

        for filepath in (filepaths[0], filepaths[-1]):
            vtk_density, vtk_velocity = self.parse_vtk(filepath, shape)
            np.testing.assert_allclose(vtk_density, node_density, rtol=5e-7)
            np.testing.assert_allclose(
                vtk_velocity, node_velocity * self.lbf.tau, rtol=5e-7)
            # TODO Walberla
            # np.testing.assert_allclose(node_pressure, vtk_pressure, rtol=5e-7)


if __name__ == '__main__':
    ut.main()
