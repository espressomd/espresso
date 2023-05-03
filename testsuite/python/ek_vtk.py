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

import espressomd
import espressomd.electrokinetics
import espressomd.shapes

with contextlib.suppress(ImportError):
    import vtk  # pylint: disable=unused-import
    import espressomd.io.vtk


class EKWalberlaWrite:
    system = espressomd.System(box_l=[6, 7, 5])
    system.time_step = 1.0
    system.cell_system.skin = 0.4

    density = 1.0
    diffusion = 0.1
    valency = 0.0
    kT = 1.0

    @classmethod
    def setUpClass(cls):
        cls.lattice = espressomd.lb.LatticeWalberla(
            n_ghost_layers=1, agrid=0.5)
        cls.solver = espressomd.electrokinetics.EKNone(lattice=cls.lattice)

    def setUp(self):
        self.species = espressomd.electrokinetics.EKSpecies(
            lattice=self.lattice, density=self.density, kT=self.kT,
            diffusion=self.diffusion, valency=self.valency,
            advection=False, friction_coupling=False, ext_efield=[0., 0., 0.],
            tau=1.0, **self.ek_params)
        self.system.ekcontainer.add(self.species)
        self.system.ekcontainer.tau = 1.0
        self.system.ekcontainer.solver = self.solver

    def tearDown(self):
        self.system.ekcontainer.clear()

    def test_vtk(self):
        '''
        Check VTK files. Keep in mind the VTK module writes in single-precision.
        '''
        dist = 1.5 * self.lattice.agrid
        self.species.add_boundary_from_shape(
            shape=espressomd.shapes.Wall(normal=[1, 0, 0], dist=dist),
            value=0.0, boundary_type=espressomd.electrokinetics.DensityBoundary)
        self.species.add_boundary_from_shape(
            shape=espressomd.shapes.Wall(
                normal=[-1, 0, 0], dist=-(self.system.box_l[0] - dist)),
            value=0.0, boundary_type=espressomd.electrokinetics.DensityBoundary)

        n_steps = 100
        ek_steps = int(np.floor(n_steps * self.system.ekcontainer.tau))
        shape = tuple(self.lattice.shape)
        shape = (shape[0] - 4, *shape[1:])
        vtk_reader = espressomd.io.vtk.VTKReader()
        label_density = 'density'

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
        ek_vtk = espressomd.electrokinetics.VTKOutput(
            species=self.species, identifier=label_vtk_continuous,
            observables=vtk_obs, delta_N=1, base_folder=str(path_vtk_root))
        ek_vtk.disable()
        ek_vtk.enable()
        self.system.integrator.run(n_steps)
        ek_vtk = espressomd.electrokinetics.VTKOutput(
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
        last_frames = []
        for filepath in (path_vtk_end, path_vtk_continuous[-1]):
            grids = vtk_reader.parse(filepath)
            last_frames.append(grids[label_density])

        # check VTK output is identical in both continuous and manual mode
        for i in range(len(last_frames[0])):
            np.testing.assert_allclose(last_frames[0][i],
                                       last_frames[1][i], atol=1e-10)

        # check VTK values match node values in the final time step
        node_density = np.copy(self.species[2:-2, :, :].density)

        for vtk_density in last_frames:
            np.testing.assert_allclose(vtk_density, node_density, rtol=5e-7)

    @ut.skipIf(system.cell_system.get_state()["n_nodes"] > 4,
               "this test is slow on more than 4 MPI ranks")
    def test_exceptions(self):
        label_invalid_obs = f'test_lb_vtk_{self.ek_vtk_id}_invalid_obs'
        error_msg = r"Only the following VTK observables are supported: \['density'\], got 'dens'"
        with self.assertRaisesRegex(ValueError, error_msg):
            espressomd.electrokinetics.VTKOutput(
                species=self.species, identifier=label_invalid_obs, delta_N=0,
                observables=['dens'])
        ek_vtk_manual_id = f'test_ek_vtk_{self.ek_vtk_id}_manual'
        ek_vtk_auto_id = f'test_ek_vtk_{self.ek_vtk_id}_auto'
        vtk_manual = espressomd.electrokinetics.VTKOutput(
            species=self.species, identifier=ek_vtk_manual_id, delta_N=0,
            observables=['density'])
        vtk_auto = espressomd.electrokinetics.VTKOutput(
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
        vtk_cleared = espressomd.electrokinetics.VTKOutput(
            species=self.species, identifier=label_cleared,
            observables=['density'])
        self.system.actors.clear()
        vtk_cleared.write()
        espressomd.electrokinetics.VTKOutput(species=self.species,
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
