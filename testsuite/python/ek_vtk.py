#
# Copyright (C) 2022 The ESPResSo project
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
import espressomd.EKSpecies
import espressomd.shapes


class EKWalberlaWrite:
    system = espressomd.System(box_l=[12, 14, 16])
    system.time_step = 1.0
    system.cell_system.skin = 0.4

    density = 1.0
    diffusion = 0.1
    valency = 0.0
    kT = 1.0

    @classmethod
    def setUpClass(cls):
        cls.lattice = espressomd.lb.LatticeWalberla(
            n_ghost_layers=1, agrid=1.0)
        cls.solver = espressomd.EKSpecies.EKNone(lattice=cls.lattice)

    def setUp(self):
        self.species = espressomd.EKSpecies.EKSpecies(
            lattice=self.lattice, density=self.density, kT=self.kT,
            diffusion=self.diffusion, valency=self.valency,
            advection=False, friction_coupling=False, ext_efield=[0., 0., 0.],
            **self.ek_params)
        self.system.ekcontainer.add(self.species)
        self.system.ekcontainer.tau = 1.0
        self.system.ekcontainer.solver = self.solver

    def tearDown(self):
        self.system.ekcontainer.clear()

    def get_cell_array(self, cell, name, shape):
        return vtk.util.numpy_support.vtk_to_numpy(
            cell.GetArray(name)).reshape(shape, order='F')

    def parse_vtk(self, filepath, shape):
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(str(filepath))
        reader.Update()

        data = reader.GetOutput()
        cell = data.GetCellData()

        vtk_density = self.get_cell_array(cell, 'density', shape)
        return vtk_density

    def test_vtk(self):
        '''
        Check VTK files. Keep in mind the VTK module writes in single-precision.
        '''
        self.species.add_boundary_from_shape(
            shape=espressomd.shapes.Wall(normal=[1, 0, 0], dist=1.5),
            value=0.0, boundary_type=espressomd.EKSpecies.DensityBoundary)
        self.species.add_boundary_from_shape(
            shape=espressomd.shapes.Wall(normal=[-1, 0, 0], dist=-10.5),
            value=0.0, boundary_type=espressomd.EKSpecies.DensityBoundary)
        shape = [8, 14, 16]
        x_offset = 2

        n_steps = 100
        ek_steps = int(np.floor(n_steps * self.system.ekcontainer.tau))

        with tempfile.TemporaryDirectory() as tmp_directory:
            path_vtk_root = pathlib.Path(tmp_directory)
            label_vtk_end = f'test_ek_vtk_{self.ek_vtk_id}_end'
            label_vtk_continuous = f'test_ek_vtk_{self.ek_vtk_id}_continuous'
            path_vtk_end = path_vtk_root / label_vtk_end / 'simulation_step_0.vtu'
            path_vtk_continuous = [
                path_vtk_root / label_vtk_continuous / f'simulation_step_{i}.vtu' for i in range(ek_steps)]
            filepaths = [path_vtk_end] + path_vtk_continuous

        # write VTK files
        vtk_obs = ['density']
        ek_vtk = espressomd.EKSpecies.EKVTKOutput(
            species=self.species, identifier=label_vtk_continuous,
            observables=vtk_obs, delta_N=1, base_folder=str(path_vtk_root))
        ek_vtk.disable()
        ek_vtk.enable()
        self.system.integrator.run(n_steps)
        ek_vtk = espressomd.EKSpecies.EKVTKOutput(
            species=self.species, identifier=label_vtk_end,
            observables=vtk_obs, delta_N=0, base_folder=str(path_vtk_root))
        ek_vtk.write()
        self.assertEqual(sorted(ek_vtk.observables), sorted(vtk_obs))
        self.assertEqual(ek_vtk.valid_observables(), {"density"})

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

        # read VTK output of final time step
        last_frame_outputs = []
        for filepath in (path_vtk_end, path_vtk_continuous[-1]):
            last_frame_outputs.append(self.parse_vtk(filepath, shape))

            # check VTK output is identical in both continuous and manual mode
        for i in range(len(last_frame_outputs[0])):
            np.testing.assert_allclose(last_frame_outputs[0][i],
                                       last_frame_outputs[1][i], atol=1e-10)

        # check VTK values match node values in the final time step
        node_density = np.copy(self.species[x_offset:-x_offset, :, :].density)

        for vtk_density in last_frame_outputs:
            np.testing.assert_allclose(vtk_density, node_density, rtol=5e-7)

    def test_exceptions(self):
        label_invalid_obs = f'test_lb_vtk_{self.ek_vtk_id}_invalid_obs'
        error_msg = r"Only the following VTK observables are supported: \['density'\], got 'dens'"
        with self.assertRaisesRegex(ValueError, error_msg):
            espressomd.EKSpecies.EKVTKOutput(
                species=self.species, identifier=label_invalid_obs, delta_N=0,
                observables=['dens'])
        ek_vtk_manual_id = f'test_ek_vtk_{self.ek_vtk_id}_manual'
        ek_vtk_auto_id = f'test_ek_vtk_{self.ek_vtk_id}_auto'
        vtk_manual = espressomd.EKSpecies.EKVTKOutput(
            species=self.species, identifier=ek_vtk_manual_id, delta_N=0,
            observables=['density'])
        vtk_auto = espressomd.EKSpecies.EKVTKOutput(
            species=self.species, identifier=ek_vtk_auto_id, delta_N=1,
            observables=['density'])
        with self.assertRaisesRegex(RuntimeError, 'Automatic VTK callbacks cannot be triggered manually'):
            vtk_auto.write()
        with self.assertRaisesRegex(RuntimeError, 'Manual VTK callbacks cannot be disabled'):
            vtk_manual.disable()
        with self.assertRaisesRegex(RuntimeError, 'Manual VTK callbacks cannot be enabled'):
            vtk_manual.enable()

        # can still use VTK when the EK actor has been cleared but not deleted
        label_cleared = f'test_ek_vtk_{self.ek_vtk_id}_cleared'
        vtk_cleared = espressomd.EKSpecies.EKVTKOutput(
            species=self.species, identifier=label_cleared,
            observables=['density'])
        self.system.actors.clear()
        vtk_cleared.write()
        espressomd.EKSpecies.EKVTKOutput(species=self.species,
                                         identifier=label_cleared + '_1',
                                         observables=['density'])


@utx.skipIfMissingModules("vtk")
@utx.skipIfMissingFeatures("WALBERLA")
class LBWalberlaWrite(EKWalberlaWrite, ut.TestCase):
    ek_params = {'single_precision': False}
    ek_vtk_id = 'double_precision'


@utx.skipIfMissingModules("vtk")
@utx.skipIfMissingFeatures("WALBERLA")
class LBWalberlaWriteSinglePrecision(EKWalberlaWrite, ut.TestCase):
    ek_params = {'single_precision': True}
    ek_vtk_id = 'single_precision'


if __name__ == '__main__':
    ut.main()
