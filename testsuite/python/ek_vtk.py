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
import espressomd.EKSpecies
import espressomd.shapes


@skipIfMissingPythonPackage
@utx.skipIfMissingFeatures("LB_WALBERLA")
class EKWalberlaWrite(ut.TestCase):
    ek_vtk_id = 'double_precision'

    system = espressomd.System(box_l=3 * [16])
    system.time_step = 1.0
    system.cell_system.skin = 0.4

    lattice = espressomd.lb.LatticeWalberla(
        box_size=system.box_l, n_ghost_layers=1, agrid=1.0)
    density = 1.0
    diffusion = 0.1
    valency = 0.0
    kT = 1.0

    solver = espressomd.EKSpecies.EKNone(lattice=lattice, permittivity=1.0)

    def setUp(self):
        self.species = espressomd.EKSpecies.EKSpecies(
            lattice=self.lattice, density=self.density, kT=self.kT,
            diffusion=self.diffusion, valency=self.valency,
            advection=False, friction_coupling=False, ext_efield=[0, 0, 0])
        self.system.ekcontainer.add(self.species)
        self.system.ekcontainer.tau = 1.0
        self.system.ekcontainer.solver = self.solver

    def tearDown(self):
        self.system.ekcontainer.clear()

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

        vtk_density = self.get_cell_array(cell, 'density', shape)
        return vtk_density

    def test_vtk(self):
        '''
        Check VTK files. Keep in mind VTK files are written with
        float precision.
        '''
        x_offset = 0
        shape = [16, 16, 16]

        self.species.add_boundary_from_shape(shape=espressomd.shapes.Wall(normal=[1, 0, 0], dist=1.5),
                                             value=0.0, boundary_type=espressomd.EKSpecies.DensityBoundary)
        self.species.add_boundary_from_shape(shape=espressomd.shapes.Wall(normal=[-1, 0, 0], dist=-14.5),
                                             value=0.0, boundary_type=espressomd.EKSpecies.DensityBoundary)
        x_offset = 2
        shape[0] = 12

        n_steps = 100
        ek_steps = int(np.floor(n_steps * self.system.ekcontainer.tau))
        label_vtk_end = f'test__vtk_{self.ek_vtk_id}_end'
        label_vtk_continuous = f'test_ek_vtk_{self.ek_vtk_id}_continuous'
        filepath_vtk_end = f'vtk_out/{label_vtk_end}/simulation_step_0.vtu'
        filepath_vtk_continuous = [
            f'vtk_out/{label_vtk_continuous}/simulation_step_{i}.vtu' for i in range(ek_steps)]
        filepaths = [filepath_vtk_end] + filepath_vtk_continuous

        # cleanup action
        self.cleanup_vtk_files(filepaths)

        # write VTK files
        vtk_obs = ['density', ]
        espressomd.EKSpecies.EKVTKOutput(species=self.species, identifier=label_vtk_continuous,
                                         observables=vtk_obs, delta_N=1)
        self.system.integrator.run(n_steps)
        ek_vtk = espressomd.EKSpecies.EKVTKOutput(species=self.species, identifier=label_vtk_end,
                                                  observables=vtk_obs, delta_N=0)
        ek_vtk.write()

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
        for i in range(shape[0]):
            for j in range(shape[1]):
                for k in range(shape[2]):
                    node = self.species[i + x_offset, j, k]
                    node_density[i, j, k] = node.density

        for vtk_density in last_frame_outputs:
            np.testing.assert_allclose(vtk_density, node_density, rtol=5e-7)

        self.cleanup_vtk_files(filepaths)

    def test_exceptions(self):
        label_invalid_obs = f'test_lb_vtk_{self.ek_vtk_id}_invalid_obs'
        error_msg = r"Only the following VTK observables are supported: \['density'\], got \['dens'\]"
        with self.assertRaisesRegex(ValueError, error_msg):
            espressomd.EKSpecies.EKVTKOutput(species=self.species, identifier=label_invalid_obs,
                                             observables=['dens'])
        label_manual_species = f'test_ek_vtk_{self.ek_vtk_id}_manual_species'
        label_auto_species = f'test_ek_vtk_{self.ek_vtk_id}_auto_species'
        vtk_manual = espressomd.EKSpecies.EKVTKOutput(species=self.species, identifier=label_manual_species,
                                                      observables=['density'], delta_N=0)
        vtk_auto = espressomd.EKSpecies.EKVTKOutput(species=self.species, identifier=label_auto_species,
                                                    observables=['density'], delta_N=1)
        with self.assertRaisesRegex(RuntimeError, 'Automatic VTK callbacks cannot be triggered manually'):
            vtk_auto.write()
        with self.assertRaisesRegex(RuntimeError, 'Manual VTK callbacks cannot be disabled'):
            vtk_manual.disable()
        with self.assertRaisesRegex(RuntimeError, 'Manual VTK callbacks cannot be enabled'):
            vtk_manual.enable()

        # can still use VTK when the LB actor has been cleared but not deleted
        label_cleared = f'test_ek_vtk_{self.ek_vtk_id}_cleared_species'
        vtk_cleared = espressomd.EKSpecies.EKVTKOutput(species=self.species, identifier=label_cleared,
                                                       observables=['density'])
        self.system.actors.clear()
        vtk_cleared.write()
        espressomd.EKSpecies.EKVTKOutput(species=self.species, observables=['density'],
                                         identifier=label_cleared + '_1')


if __name__ == '__main__':
    ut.main()
