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
import espressomd.shapes


class LBWrite:
    system = espressomd.System(box_l=3 * [16])
    system.time_step = 0.01
    system.cell_system.skin = 0.4

    def setUp(self):
        self.lbf = self.lb_class(
            kT=0, agrid=1.0, density=1.0, viscosity=1.0, tau=0.1,
            ext_force_density=[0, 0.03, 0], **self.lb_params)
        self.system.actors.add(self.lbf)

    def tearDown(self):
        if self.lbf is not None:
            self.lbf.clear_boundaries()
        self.system.actors.clear()
        self.system.thermostat.turn_off()

    def get_vtk_folder_names(self, filepaths):
        return set(os.path.dirname(filepath) for filepath in filepaths)

    def cleanup_vtk_files(self, filepaths):
        '''
        Remove VTK files (.vtk and .vtu files), their folder (if empty) and
        their summaries (.vtd files).
        '''
        for filepath in filepaths:
            if os.path.exists(filepath):
                os.remove(filepath)
        for dirname in self.get_vtk_folder_names(filepaths):
            if os.path.exists(dirname) and len(os.listdir(dirname)) == 0:
                os.rmdir(dirname)
            filepath = dirname + '.pvd'
            if os.path.exists(filepath):
                os.remove(filepath)

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
            cell, 'VelocityFromVelocityField', shape + [3])
        vtk_pressure = self.get_cell_array(
            cell, 'PressureTensorFromPDF', shape + [3, 3])
        return vtk_density, vtk_velocity, vtk_pressure

    def test_vtk(self):
        '''
        Check VTK files. Keep in mind VTK files are written with
        float precision.
        '''
        x_offset = 0
        shape = [16, 16, 16]
        self.lbf.add_boundary_from_shape(
            espressomd.shapes.Wall(normal=[1, 0, 0], dist=1.5))
        self.lbf.add_boundary_from_shape(
            espressomd.shapes.Wall(normal=[-1, 0, 0], dist=-14.5))
        x_offset = 2
        shape[0] = 12

        n_steps = 100
        lb_steps = int(np.floor(n_steps * self.lbf.tau))
        label_vtk_end = f'test_lb_vtk_{self.lb_vtk_id}_end'
        label_vtk_continuous = f'test_lb_vtk_{self.lb_vtk_id}_continuous'
        filepath_vtk_end = f'vtk_out/{label_vtk_end}/simulation_step_0.vtu'
        filepath_vtk_continuous = [
            f'vtk_out/{label_vtk_continuous}/simulation_step_{i}.vtu' for i in range(lb_steps)]
        filepaths = [filepath_vtk_end] + filepath_vtk_continuous

        # cleanup action
        self.cleanup_vtk_files(filepaths)

        # write VTK files
        vtk_obs = ['density', 'velocity_vector', 'pressure_tensor']
        self.lbf.add_vtk_writer(label_vtk_continuous, vtk_obs, delta_N=1)
        self.system.integrator.run(n_steps)
        lb_vtk = self.lbf.add_vtk_writer(label_vtk_end, vtk_obs, delta_N=0)
        lb_vtk.write()

        # check VTK files exist
        for filepath in filepaths:
            self.assertTrue(
                os.path.exists(filepath),
                f'VTK file "{filepath}" not written to disk')
        for dirname in self.get_vtk_folder_names(filepaths):
            filepath = dirname + '.pvd'
            self.assertTrue(
                os.path.exists(filepath),
                f'VTK summary file "{filepath}" not written to disk')

        # check velocity profile is symmetric at all time steps
        for filepath in filepaths:
            vtk_velocity = self.parse_vtk(filepath, shape)[1]
            v_profile = np.mean(
                np.linalg.norm(vtk_velocity, axis=-1),
                axis=(1, 2))
            np.testing.assert_allclose(v_profile, v_profile[::-1], atol=5e-6)

        # check pressure tensor is symmetric at all time steps
        for filepath in filepaths:
            vtk_velocity = self.parse_vtk(filepath, shape)[1]
            v_profile = np.mean(
                np.linalg.norm(vtk_velocity, axis=-1),
                axis=(1, 2))
            np.testing.assert_allclose(v_profile, v_profile[::-1], atol=5e-6)

        # read VTK output of final time step
        last_frame_outputs = []
        for filepath in (filepath_vtk_end, filepath_vtk_continuous[-1]):
            last_frame_outputs.append(self.parse_vtk(filepath, shape))

        # check the VTK output is identical in both continuous and manual mode
        for i in range(len(last_frame_outputs[0])):
            np.testing.assert_allclose(last_frame_outputs[0][i],
                                       last_frame_outputs[1][i], atol=1e-10)

        # check VTK values match node values in the final time step
        node_density = np.zeros(shape)
        node_velocity = np.zeros(shape + [3])
        node_pressure = np.zeros(shape + [3, 3])
        tau = self.lbf.tau
        for i in range(shape[0]):
            for j in range(shape[1]):
                for k in range(shape[2]):
                    node = self.lbf[i + x_offset, j, k]
                    node_density[i, j, k] = node.density
                    node_velocity[i, j, k] = node.velocity * tau
                    node_pressure[i, j, k] = node.pressure_tensor * tau**2

        for vtk_density, vtk_velocity, vtk_pressure in last_frame_outputs:
            np.testing.assert_allclose(vtk_density, node_density, rtol=5e-7)
            np.testing.assert_allclose(vtk_velocity, node_velocity, rtol=5e-7)
            # TODO WALBERLA mismatch in off-diagonal terms
            np.testing.assert_allclose(vtk_pressure, node_pressure, atol=1e-3)

        self.cleanup_vtk_files(filepaths)

    def test_exceptions(self):
        label_invalid_obs = f'test_lb_vtk_{self.lb_vtk_id}_invalid_obs'
        error_msg = r"Only the following VTK observables are supported: \['density', 'pressure_tensor', 'velocity_vector'\], got \['dens'\]"
        with self.assertRaisesRegex(ValueError, error_msg):
            self.lbf.add_vtk_writer(label_invalid_obs, ['dens'])
        label_manual_lbf = f'test_lb_vtk_{self.lb_vtk_id}_manual_lbf'
        label_auto_lbf = f'test_lb_vtk_{self.lb_vtk_id}_auto_lbf'
        vtk_manual = self.lbf.add_vtk_writer(label_manual_lbf, ['density'])
        vtk_auto = self.lbf.add_vtk_writer(
            label_auto_lbf, ['density'], delta_N=1)
        with self.assertRaisesRegex(RuntimeError, 'Automatic VTK callbacks cannot be triggered manually'):
            vtk_auto.write()
        with self.assertRaisesRegex(RuntimeError, 'Manual VTK callbacks cannot be disabled'):
            vtk_manual.disable()
        with self.assertRaisesRegex(RuntimeError, 'Manual VTK callbacks cannot be enabled'):
            vtk_manual.enable()

        # can still use VTK when the LB actor has been cleared but not deleted
        label_cleared = f'test_lb_vtk_{self.lb_vtk_id}_cleared_lbf'
        vtk_cleared = self.lbf.add_vtk_writer(label_cleared, ['density'])
        self.system.actors.clear()
        vtk_cleared.write()
        espressomd.lb.VTKOutput(lb_fluid=self.lbf, observables=['density'],
                                identifier=label_cleared + '_1')

        # cannot use VTK when the LB actor has expired
        label_expired = f'test_lb_vtk_{self.lb_vtk_id}_expired_lbf'
        vtk_expired = self.lbf.add_vtk_writer(label_expired, ['density'])
        self.lbf = None
        with self.assertRaisesRegex(RuntimeError, 'Attempted access to uninitialized LBWalberla instance'):
            vtk_expired.write()


@skipIfMissingPythonPackage
@utx.skipIfMissingFeatures("LB_WALBERLA")
class LBWalberlaWrite(LBWrite, ut.TestCase):
    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': False}
    lb_vtk_id = 'double_precision'


@skipIfMissingPythonPackage
@utx.skipIfMissingFeatures("LB_WALBERLA")
class LBWalberlaWriteSinglePrecision(LBWrite, ut.TestCase):
    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': True}
    lb_vtk_id = 'single_precision'


if __name__ == '__main__':
    ut.main()
