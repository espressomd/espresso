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

import pathlib
import tempfile
import contextlib
import numpy as np

with contextlib.suppress(ImportError):
    import vtk
    import vtk.util.numpy_support

import espressomd
import espressomd.lb
import espressomd.shapes


class TestLBWrite:
    """
    Set up a planar Poiseuille flow and write fluid to VTK files.
    """

    system = espressomd.System(box_l=[12, 14, 16])
    system.time_step = 0.01
    system.cell_system.skin = 0.4

    def setUp(self):
        self.lbf = self.lb_class(
            kT=0, agrid=1.0, density=1.0, viscosity=1.0, tau=0.1,
            ext_force_density=[0, 0.03, 0], **self.lb_params)
        self.system.actors.add(self.lbf)

    def tearDown(self):
        self.system.actors.clear()
        self.system.thermostat.turn_off()

    def get_cell_array(self, cell, name, shape):
        return vtk.util.numpy_support.vtk_to_numpy(
            cell.GetArray(name)).reshape(shape, order='F')

    def parse_vtk(self, filepath, shape):
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(str(filepath))
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
        Check VTK files. Keep in mind the VTK module writes in single-precision.
        '''
        self.lbf.add_boundary_from_shape(
            espressomd.shapes.Wall(normal=[1, 0, 0], dist=1.5))
        self.lbf.add_boundary_from_shape(
            espressomd.shapes.Wall(normal=[-1, 0, 0], dist=-10.5))
        shape = [8, 14, 16]

        n_steps = 100
        lb_steps = int(np.floor(n_steps * self.lbf.tau))

        with tempfile.TemporaryDirectory() as tmp_directory:
            path_vtk_root = pathlib.Path(tmp_directory)
            label_vtk_end = f'test_lb_vtk_{self.lb_vtk_id}_end'
            label_vtk_continuous = f'test_lb_vtk_{self.lb_vtk_id}_continuous'
            path_vtk_end = path_vtk_root / label_vtk_end / 'simulation_step_0.vtu'
            path_vtk_continuous = [
                path_vtk_root / label_vtk_continuous / f'simulation_step_{i}.vtu' for i in range(lb_steps)]
            filepaths = [path_vtk_end] + path_vtk_continuous

            # write VTK files
            vtk_obs = ['density', 'velocity_vector', 'pressure_tensor']
            lb_vtk = espressomd.lb.VTKOutput(
                lb_fluid=self.lbf, identifier=label_vtk_continuous, delta_N=1,
                observables=vtk_obs, base_folder=str(path_vtk_root))
            lb_vtk.disable()
            lb_vtk.enable()
            self.system.integrator.run(n_steps)
            lb_vtk = espressomd.lb.VTKOutput(
                lb_fluid=self.lbf, identifier=label_vtk_end, delta_N=0,
                observables=vtk_obs, base_folder=str(path_vtk_root))
            lb_vtk.write()
            self.assertEqual(sorted(lb_vtk.observables), sorted(vtk_obs))
            self.assertEqual(lb_vtk.valid_observables(),
                             {"density", "pressure_tensor", "velocity_vector"})

            # check VTK files exist
            for filepath in filepaths:
                self.assertTrue(
                    filepath.exists(),
                    f'VTK file "{filepath}" not written to disk')
            for filepath in [path_vtk_end.parent.with_suffix('.pvd'),
                             path_vtk_continuous[0].parent.with_suffix('.pvd')]:
                self.assertTrue(
                    filepath.exists(),
                    f'VTK summary file "{filepath}" not written to disk')

            # check velocity profile is symmetric at all time steps
            for filepath in filepaths:
                vtk_velocity = self.parse_vtk(filepath, shape)[1]
                v_profile = np.mean(
                    np.linalg.norm(vtk_velocity, axis=-1),
                    axis=(1, 2))
                np.testing.assert_allclose(
                    v_profile, v_profile[::-1], atol=1e-7)

            # check scalar pressure is symmetric at all time steps
            for filepath in filepaths:
                vtk_pressure = self.parse_vtk(filepath, shape)[2]
                p_profile = np.mean(
                    np.trace(vtk_pressure, axis1=-2, axis2=-1),
                    axis=(1, 2))
                np.testing.assert_allclose(
                    p_profile, p_profile[::-1], atol=1e-6)

            # read VTK output of final time step
            last_frame_outputs = []
            for filepath in (path_vtk_end, path_vtk_continuous[-1]):
                last_frame_outputs.append(self.parse_vtk(filepath, shape))

            # check VTK output is identical in both continuous and manual mode
            for i in range(len(last_frame_outputs[0])):
                np.testing.assert_allclose(last_frame_outputs[0][i],
                                           last_frame_outputs[1][i], atol=1e-10)

            # check VTK values match node values in the final time step
            tau = self.lbf.tau
            lb_density = np.copy(self.lbf[2:-2, :, :].density)
            lb_velocity = np.copy(self.lbf[2:-2, :, :].velocity) * tau
            lb_pressure = np.copy(
                self.lbf[2:-2, :, :].pressure_tensor) * tau**2

            for vtk_density, vtk_velocity, vtk_pressure in last_frame_outputs:
                np.testing.assert_allclose(
                    vtk_density, lb_density, rtol=1e-10, atol=0.)
                np.testing.assert_allclose(
                    vtk_velocity, lb_velocity, rtol=1e-7, atol=0.)
                np.testing.assert_allclose(
                    vtk_pressure, lb_pressure, rtol=1e-7, atol=0.)

    def test_exceptions(self):
        label_invalid_obs = f'test_lb_vtk_{self.lb_vtk_id}_invalid_obs'
        error_msg = r"Only the following VTK observables are supported: \['density', 'pressure_tensor', 'velocity_vector'\], got 'dens'"
        with self.assertRaisesRegex(ValueError, error_msg):
            espressomd.lb.VTKOutput(
                lb_fluid=self.lbf, identifier=label_invalid_obs, delta_N=0,
                observables=['dens'])
        lb_vtk_manual_id = f'test_lb_vtk_{self.lb_vtk_id}_manual'
        lb_vtk_auto_id = f'test_lb_vtk_{self.lb_vtk_id}_auto'
        vtk_manual = espressomd.lb.VTKOutput(
            lb_fluid=self.lbf, identifier=lb_vtk_manual_id, delta_N=0,
            observables=['density'])
        vtk_auto = espressomd.lb.VTKOutput(
            lb_fluid=self.lbf, identifier=lb_vtk_auto_id, delta_N=1,
            observables=['density'])
        with self.assertRaisesRegex(RuntimeError, 'Automatic VTK callbacks cannot be triggered manually'):
            vtk_auto.write()
        with self.assertRaisesRegex(RuntimeError, 'Manual VTK callbacks cannot be disabled'):
            vtk_manual.disable()
        with self.assertRaisesRegex(RuntimeError, 'Manual VTK callbacks cannot be enabled'):
            vtk_manual.enable()

        # can still use VTK when the LB actor has been cleared but not deleted
        label_cleared = f'test_lb_vtk_{self.lb_vtk_id}_cleared'
        vtk_cleared = espressomd.lb.VTKOutput(
            lb_fluid=self.lbf, identifier=label_cleared,
            observables=['density'])
        self.system.actors.clear()
        vtk_cleared.write()
        espressomd.lb.VTKOutput(lb_fluid=self.lbf,
                                identifier=label_cleared + '_1',
                                observables=['density'])

        # cannot use VTK when the LB actor has expired
        label_expired = f'test_lb_vtk_{self.lb_vtk_id}_expired_lbf'
        vtk_expired = espressomd.lb.VTKOutput(
            lb_fluid=self.lbf, identifier=label_expired, observables=['density'])
        self.lbf = None
        with self.assertRaisesRegex(RuntimeError, 'Attempted access to uninitialized LBWalberla instance'):
            vtk_expired.write()

        # cannot use VTK when there are no LB objects available
        label_unavailable = f'test_lb_vtk_{self.lb_vtk_id}_unavailable_lbf'

        class VTKOutputWithoutLbfluid(espressomd.lb.VTKOutput):
            def validate_params(self, params):
                pass
        with self.assertRaisesRegex(RuntimeError, 'Attempted access to uninitialized LBWalberla instance'):
            VTKOutputWithoutLbfluid(identifier=label_unavailable,
                                    observables=['density'])


@utx.skipIfMissingModules("vtk")
@utx.skipIfMissingFeatures("WALBERLA")
class LBWalberlaWrite(TestLBWrite, ut.TestCase):
    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': False}
    lb_vtk_id = 'double_precision'


@utx.skipIfMissingModules("vtk")
@utx.skipIfMissingFeatures("WALBERLA")
class LBWalberlaWriteSinglePrecision(TestLBWrite, ut.TestCase):
    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {'single_precision': True}
    lb_vtk_id = 'single_precision'


if __name__ == '__main__':
    ut.main()
