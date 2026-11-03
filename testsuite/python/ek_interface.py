#
# Copyright (C) 2010-2026 The ESPResSo project
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

import numpy as np
import itertools
import pathlib
import tempfile
import unittest as ut
import unittest_decorators as utx

import espressomd
import espressomd.lb
import espressomd.electrokinetics


class EKTest:

    """
    Basic tests for the electrokinetics implementation

    """
    system = espressomd.System(box_l=3 * [6.0])
    np.random.seed(1)
    params = {"tau": 0.01, "agrid": 0.5}
    ek_species_params = {"kT": 1.5,
                         "density": 0.85,
                         "valency": 0.0,
                         "diffusion": 0.1,
                         "advection": False,
                         "friction_coupling": False,
                         "ext_efield": [0.1, 0.2, 0.3],
                         "tau": params["tau"]}

    system.periodicity = [True, True, True]
    system.time_step = params["tau"]
    system.cell_system.skin = 1.0

    def setUp(self):
        self.system.box_l = 3 * [6.0]
        self.lattice = self.ek_lattice_class(
            n_ghost_layers=2, agrid=self.params["agrid"])
        ek_solver = espressomd.electrokinetics.EKNone(
            lattice=self.lattice, tau=self.params["tau"], **self.ek_params)
        self.system.ekcontainer = espressomd.electrokinetics.EKContainer(
            tau=self.system.time_step, solver=ek_solver)

    def tearDown(self):
        self.system.lb = None
        self.system.ekcontainer = None
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.time_step = self.params["tau"]

    def make_default_ek_species(self, **kwargs):
        return self.ek_species_class(
            lattice=self.lattice,
            **(self.ek_params | self.ek_species_params | kwargs))

    def test_ek_species(self):
        # inactive species
        ek_species = self.make_default_ek_species()
        self.check_ek_species_properties(ek_species)

        self.assertAlmostEqual(
            self.system.ekcontainer.tau,
            self.system.time_step,
            delta=self.atol)

        # activated species
        ek_species = self.make_default_ek_species()
        self.system.ekcontainer.add(ek_species)
        self.check_ek_species_properties(ek_species)

        # deactivated species
        ek_species = self.make_default_ek_species()
        self.system.ekcontainer.add(ek_species)
        self.system.ekcontainer.remove(ek_species)
        self.check_ek_species_properties(ek_species)

        # reactive species
        ek_species = self.make_default_ek_species()
        espressomd.electrokinetics.EKReactant(
            ekspecies=ek_species, stoech_coeff=-2.0, order=2.0)
        self.check_ek_species_properties(ek_species)

        # thermalized species
        ek_species = self.make_default_ek_species(thermalized=True, seed=42)
        self.assertTrue(ek_species.thermalized)
        self.assertEqual(ek_species.seed, 42)
        self.assertEqual(ek_species.rng_state, 0)
        ek_species.rng_state = 5
        self.assertEqual(ek_species.rng_state, 5)

    def check_ek_species_properties(self, species):
        agrid = self.params["agrid"]
        # check getters
        self.assertEqual(species.lattice.n_ghost_layers, 2)
        self.assertAlmostEqual(species.lattice.agrid, agrid, delta=self.atol)
        self.assertAlmostEqual(species.diffusion, 0.1, delta=self.atol)
        self.assertAlmostEqual(species.valency, 0.0, delta=self.atol)
        self.assertAlmostEqual(species.kT, 1.5, delta=self.atol)
        np.testing.assert_allclose(
            np.copy(species.ext_efield), [0.1, 0.2, 0.3], atol=self.atol)
        self.assertFalse(species.advection)
        self.assertFalse(species.friction_coupling)
        self.assertEqual(
            species.single_precision,
            self.ek_params["single_precision"])
        # check setters
        species.diffusion = 0.2
        species.valency = 0.3
        species.kT = 0.4
        ext_f = [0.01, 0.02, 0.03]
        species.ext_efield = ext_f
        species.advection = True
        species.friction_coupling = True
        self.assertAlmostEqual(species.diffusion, 0.2, delta=self.atol)
        self.assertAlmostEqual(species.valency, 0.3, delta=self.atol)
        self.assertAlmostEqual(species.kT, 0.4, delta=self.atol)
        np.testing.assert_allclose(
            np.copy(species.ext_efield), ext_f, atol=self.atol)
        self.assertTrue(species.advection)
        self.assertTrue(species.friction_coupling)
        # check node getters/setters
        self.assertAlmostEqual(species[0, 0, 0].density, 0.85, delta=self.atol)
        species[0, 0, 0].density = 0.90
        self.assertAlmostEqual(species[0, 0, 0].density, 0.90, delta=self.atol)
        with self.assertRaises(RuntimeError):
            species[0, 0, 0].density = [1, 2]
        with self.assertRaises(TypeError):
            species[0, 1].density = 1.
        # check boundary conditions
        node = species[1, 1, 1]
        self.assertIsNone(node.density_boundary)
        self.assertIsNone(node.flux_boundary)
        self.assertFalse(node.is_boundary)
        node.flux_boundary = espressomd.electrokinetics.FluxBoundary(
            [1., 2., 3.])
        self.assertIsInstance(
            node.flux_boundary,
            espressomd.electrokinetics.FluxBoundary)
        np.testing.assert_allclose(
            np.copy(node.flux_boundary.flux), [1., 2., 3.], atol=self.atol)
        self.assertTrue(node.is_boundary)
        node.density_boundary = espressomd.electrokinetics.DensityBoundary(4.5)
        self.assertIsInstance(
            node.density_boundary,
            espressomd.electrokinetics.DensityBoundary)
        np.testing.assert_allclose(
            np.copy(node.density_boundary.density), 4.5, atol=self.atol)
        self.assertTrue(node.is_boundary)
        node.density_boundary = None
        self.assertTrue(node.is_boundary)
        node.flux_boundary = None
        self.assertFalse(node.is_boundary)
        self.assertIsNone(node.density_boundary)
        self.assertIsNone(node.flux_boundary)
        with self.assertRaisesRegex(TypeError, "must be an instance of DensityBoundary or None"):
            node.density_boundary = 4.6
        with self.assertRaisesRegex(TypeError, "must be an instance of FluxBoundary or None"):
            node.flux_boundary = 4.6

    @utx.skipIfMissingFeatures(["WALBERLA_FFT"])
    def test_ek_fft_solver(self):
        ek_solver = self.ek_solver_class(
            lattice=self.lattice, permittivity=0.01,
            tau=self.params["tau"], **self.ek_params)
        self.assertEqual(ek_solver.lattice, self.lattice)
        self.assertEqual(
            ek_solver.single_precision,
            self.ek_params["single_precision"])
        self.assertAlmostEqual(ek_solver.permittivity, 0.01, delta=self.atol)
        ek_solver.permittivity = 0.05
        self.assertAlmostEqual(ek_solver.permittivity, 0.05, delta=self.atol)
        self.assertAlmostEqual(
            ek_solver.tau, self.params["tau"], delta=self.atol)

        self.system.ekcontainer.solver = ek_solver
        self.assertIsInstance(self.system.ekcontainer.solver,
                              self.ek_solver_class)
        self.assertEqual(self.system.ekcontainer.solver, ek_solver)

        ek_node = ek_solver[2, 2, 2]
        np.testing.assert_array_equal(np.copy(ek_node._index), [2, 2, 2])
        ek_node.call_method("override_index", index=[0, 1, 2])
        np.testing.assert_array_equal(np.copy(ek_node._index), [0, 1, 2])

        ek_slice = ek_solver[0:2, -4:-1, 1:2]
        np.testing.assert_array_equal(
            np.copy(ek_slice.call_method("get_slice_size")), [2, 3, 1])
        np.testing.assert_array_equal(
            np.copy(ek_slice.call_method("get_slice_ranges")),
            [[0, 8, 1], [2, 11, 2]])
        np.testing.assert_array_equal(
            np.copy(ek_slice.call_method("get_value_shape", name="potential")),
            [1])
        self.assertEqual(ek_slice.call_method("get_ek_solver_sip"), ek_solver)
        with self.assertRaisesRegex(ValueError, "Unknown Poisson solver property 'unknown'"):
            ek_slice.call_method("get_value_shape", name="unknown")
        with self.assertRaisesRegex(
                RuntimeError, "Setting potential is not supported by EKFFT"):
            ek_node.potential = 0.1
        with self.assertRaisesRegex(
                RuntimeError, "Setting potential is not supported by EKFFT"):
            ek_slice.potential = 0.1

    def test_ek_none_solver(self):
        is_gpu = self.ek_params["gpu"]
        self.assertEqual(self.system.ekcontainer.gpu, is_gpu)
        self.assertEqual(self.system.ekcontainer.solver.gpu, is_gpu)
        ek_solver = espressomd.electrokinetics.EKNone(
            lattice=self.lattice, tau=self.params["tau"], **self.ek_params)
        self.assertEqual(
            ek_solver.single_precision,
            self.ek_params["single_precision"])
        self.assertAlmostEqual(
            ek_solver.tau, self.params["tau"], delta=self.atol)

        self.system.ekcontainer.solver = ek_solver
        self.assertIsInstance(self.system.ekcontainer.solver,
                              espressomd.electrokinetics.EKNone)
        self.assertEqual(self.system.ekcontainer.solver, ek_solver)
        self.assertIsNone(ek_solver.call_method("unknown"))

        np.testing.assert_allclose(
            np.copy(ek_solver[:, :, :].potential), 0., atol=self.atol)
        self.assertAlmostEqual(
            np.copy(ek_solver[0, 0, 0].potential), 0., delta=self.atol)

        ek_solver[0, 0, 0].potential = 0.123
        self.assertAlmostEqual(
            np.copy(ek_solver[0, 0, 0].potential), 0.123, delta=self.atol)
        ek_solver[1:3, 2:4, 3:5].potential = -0.5
        np.testing.assert_allclose(
            np.copy(ek_solver[1:3, 2:4, 3:5].potential), -0.5, atol=self.atol)

        with tempfile.TemporaryDirectory() as tmpdir:
            checkpoint = pathlib.Path(tmpdir) / "eknone.cpt"
            ek_solver.save_checkpoint(path=checkpoint, binary=True)
            ek_solver[:, :, :].potential = 0.0
            np.testing.assert_allclose(
                np.copy(ek_solver[:, :, :].potential), 0., atol=self.atol)
            ek_solver.load_checkpoint(path=checkpoint, binary=True)
            self.assertAlmostEqual(
                np.copy(ek_solver[0, 0, 0].potential), 0.123, delta=self.atol)
            np.testing.assert_allclose(
                np.copy(ek_solver[1:3, 2:4, 3:5].potential), -0.5,
                atol=self.atol)

    def test_ek_species_exceptions(self):
        ek_species = self.make_default_ek_species()
        with self.assertRaisesRegex(ValueError, "Parameter 'kT' must be >= 0"):
            ek_species.kT = -0.4
        with self.assertRaisesRegex(ValueError, "Parameter 'rng_state' must be >= 0"):
            ek_species.rng_state = -2
        with self.assertRaisesRegex(RuntimeError, "This EK instance is unthermalized"):
            ek_species.rng_state = 5
        incompatible_lattice = self.ek_lattice_class(
            n_ghost_layers=1, agrid=self.params["agrid"],
            blocks_per_mpi_rank=[2, 1, 1])
        with self.assertRaisesRegex(NotImplementedError, "Using more than one block per MPI rank is not supported for EKSpecies"):
            self.ek_species_class(
                lattice=incompatible_lattice,
                **self.ek_params,
                **self.ek_species_params)
        incompatible_lattice = self.ek_lattice_class(
            n_ghost_layers=1, agrid=self.params["agrid"],
            blocks_per_mpi_rank=[1, 1, 1])
        ek_small_gl_species = self.ek_species_class(
            lattice=incompatible_lattice,
            **self.ek_params,
            **self.ek_species_params)
        if np.max(self.system.cell_system.node_grid) > 1:
            with self.assertRaisesRegex(RuntimeError, "The number of ghostlayers should be > 1 when using flux boundaries and MPI"):
                ek_small_gl_species[0, 0, 0].flux_boundary = espressomd.electrokinetics.FluxBoundary([
                    1., 2., 3.])
            with self.assertRaisesRegex(RuntimeError, "The number of ghostlayers should be > 1 when using flux boundaries and MPI"):
                ek_small_gl_species[:, :, 0].flux_boundary = espressomd.electrokinetics.FluxBoundary([
                    1., 2., 3.])
            wall_shape = espressomd.shapes.Wall(normal=[1., 0., 0.], dist=2.5)
            with self.assertRaisesRegex(RuntimeError, "The number of ghostlayers should be > 1 when using flux boundaries and MPI"):
                ek_small_gl_species.add_boundary_from_shape(
                    shape=wall_shape, value=[1., 2., 3.], boundary_type=espressomd.electrokinetics.FluxBoundary)

        self.system.ekcontainer.add(ek_species)
        with self.assertRaisesRegex(RuntimeError, "This object is already present in the list"):
            self.system.ekcontainer.add(ek_species)
        self.system.ekcontainer.remove(ek_species)
        with self.assertRaisesRegex(RuntimeError, "This object is absent from the list"):
            self.system.ekcontainer.remove(ek_species)

    @utx.skipIfMissingGPU()
    def test_ek_container_exceptions(self):
        # cannot add a CPU species if the solver is on the GPU
        ek_species_incompatible = self.make_default_ek_species(
            gpu=not self.ek_params["gpu"])
        with self.assertRaisesRegex(RuntimeError, "The EK species and the EK solver need to all be on either the CPU or GPU"):
            self.system.ekcontainer.add(ek_species_incompatible)
        self.assertEqual(len(self.system.ekcontainer), 0)
        # cannot replace a CPU solver by a GPU solver if there is at least
        # one species in the container
        self.system.ekcontainer.add(self.make_default_ek_species())
        ek_solver_incompatible = espressomd.electrokinetics.EKNone(
            lattice=self.lattice, tau=self.params["tau"],
            single_precision=self.ek_params["single_precision"],
            gpu=not self.ek_params["gpu"])
        old_ek_solver = self.system.ekcontainer.solver
        with self.assertRaisesRegex(RuntimeError, "The EK species and the EK solver need to all be on either the CPU or GPU"):
            self.system.ekcontainer.solver = ek_solver_incompatible
        self.assertEqual(self.system.ekcontainer.solver, old_ek_solver)
        # if no species are present, a GPU container can become a CPU one
        self.system.ekcontainer.clear()
        self.system.ekcontainer.solver = ek_solver_incompatible
        self.system.ekcontainer.add(ek_species_incompatible)
        self.assertEqual(
            self.system.ekcontainer.solver,
            ek_solver_incompatible)
        self.assertEqual(len(self.system.ekcontainer), 1)
        # the tau value must match for units conversion
        with self.assertRaisesRegex(RuntimeError, "EK species and EK container must have the same tau"):
            self.system.ekcontainer.add(self.make_default_ek_species(tau=200.))
        self.assertEqual(len(self.system.ekcontainer), 1)

    def test_ek_solver_exceptions(self):
        ek_solver = self.system.ekcontainer.solver
        ek_species = self.make_default_ek_species()
        self.system.ekcontainer.add(ek_species)
        incompatible_lattice = self.ek_lattice_class(
            n_ghost_layers=3, agrid=self.params["agrid"])
        incompatible_ek_solver = espressomd.electrokinetics.EKNone(
            lattice=incompatible_lattice, tau=self.params["tau"],
            **self.ek_params)
        incompatible_ek_species = self.ek_species_class(
            lattice=incompatible_lattice, **self.ek_params,
            **self.ek_species_params)
        with self.assertRaisesRegex(RuntimeError, "EKSpecies lattice incompatible with existing Poisson solver lattice"):
            self.system.ekcontainer.add(incompatible_ek_species)
        self.assertEqual(len(self.system.ekcontainer), 1)
        self.system.ekcontainer.clear()
        with self.assertRaisesRegex(RuntimeError, "Poisson solver lattice incompatible with existing EKSpecies lattice"):
            self.system.ekcontainer.solver = incompatible_ek_solver
            self.system.ekcontainer.add(incompatible_ek_species)
            self.system.ekcontainer.solver = ek_solver
        self.assertEqual(
            self.system.ekcontainer.solver, incompatible_ek_solver)
        incompatible_lattice = self.ek_lattice_class(
            n_ghost_layers=2, agrid=self.params["agrid"],
            blocks_per_mpi_rank=[2, 1, 1])
        with self.assertRaisesRegex(NotImplementedError, "Using more than one block per MPI rank is not supported for EKNone"):
            espressomd.electrokinetics.EKNone(
                lattice=incompatible_lattice, tau=self.params["tau"],
                **self.ek_params)
        with self.assertRaisesRegex(ValueError, "EK solver is of the wrong type"):
            self.system.ekcontainer.solver = ek_species
        with self.assertRaisesRegex(RuntimeError, "Parameter 'solver' is required; use EKNone if all species are electrically neutral"):
            espressomd.electrokinetics.EKContainer(tau=1.)

        if espressomd.has_features("WALBERLA_FFT"):
            ek_solver = self.ek_solver_class(
                lattice=self.lattice, permittivity=0.01, tau=self.params["tau"], **self.ek_params)
            with self.assertRaisesRegex(NotImplementedError, "Cannot serialize EK Poisson solver node objects"):
                ek_solver[0, 0, 0].__reduce__()
            with self.assertRaisesRegex(NotImplementedError, "Cannot serialize EK Poisson solver slice objects"):
                ek_solver[0:1, 0:1, 0:1].__reduce__()

            solver_sp = self.ek_solver_class(
                lattice=self.lattice, permittivity=0.1, tau=self.params["tau"], single_precision=True)
            solver_dp = self.ek_solver_class(
                lattice=self.lattice, permittivity=0.1, tau=self.params["tau"], single_precision=False)
            species_sp = self.make_default_ek_species(single_precision=True)
            species_dp = self.make_default_ek_species(single_precision=False)
            self.system.ekcontainer.clear()
            # EKNone has no effect and its floating-point precision is ignored
            solver_none_sp = espressomd.electrokinetics.EKNone(
                lattice=self.lattice, tau=self.params["tau"], single_precision=True,
                gpu=self.ek_params["gpu"])
            solver_none_dp = espressomd.electrokinetics.EKNone(
                lattice=self.lattice, tau=self.params["tau"], single_precision=True,
                gpu=self.ek_params["gpu"])
            self.system.ekcontainer.solver = solver_none_sp
            self.system.ekcontainer.add(species_dp)  # mismatch allowed
            self.system.ekcontainer.clear()
            self.system.ekcontainer.solver = solver_none_dp
            self.system.ekcontainer.add(species_sp)  # mismatch allowed
            self.system.ekcontainer.clear()
            # FFT solvers check floating-point precision of species
            self.system.ekcontainer.solver = solver_sp
            with self.assertRaisesRegex(RuntimeError, "Cannot mix single and double precision kernels"):
                self.system.ekcontainer.add(species_dp)
            self.assertEqual(len(self.system.ekcontainer), 0)
            self.system.ekcontainer.solver = solver_dp
            with self.assertRaisesRegex(RuntimeError, "Cannot mix single and double precision kernels"):
                self.system.ekcontainer.add(species_sp)
            self.assertEqual(len(self.system.ekcontainer), 0)
            # EKSpecies check floating-point precision of new species/solvers
            self.system.ekcontainer.solver = solver_none_sp
            self.system.ekcontainer.add(species_dp)
            with self.assertRaisesRegex(RuntimeError, "Cannot mix single and double precision kernels"):
                self.system.ekcontainer.add(species_sp)
            self.assertEqual(len(self.system.ekcontainer), 1)
            with self.assertRaisesRegex(RuntimeError, "Cannot mix single and double precision kernels"):
                self.system.ekcontainer.solver = solver_sp
            self.assertEqual(self.system.ekcontainer.solver, solver_none_sp)
            self.system.ekcontainer.clear()
            self.system.ekcontainer.add(species_sp)
            with self.assertRaisesRegex(RuntimeError, "Cannot mix single and double precision kernels"):
                self.system.ekcontainer.add(species_dp)
            self.assertEqual(len(self.system.ekcontainer), 1)
            with self.assertRaisesRegex(RuntimeError, "Cannot mix single and double precision kernels"):
                self.system.ekcontainer.solver = solver_dp
            self.assertEqual(self.system.ekcontainer.solver, solver_none_sp)
            self.system.ekcontainer.clear()

    def test_parameter_change_exceptions(self):
        ek_solver = self.system.ekcontainer.solver
        ek_species = self.make_default_ek_species()
        self.system.ekcontainer.add(ek_species)
        self.system.ekcontainer.solver = ek_solver
        with self.assertRaisesRegex(Exception, "Temperature change not supported by EK"):
            self.system.thermostat.set_langevin(kT=1., seed=42, gamma=1.)
        with self.assertRaisesRegex(ValueError, "must be an integer multiple of the MD time_step"):
            self.system.time_step /= 1.7
        self.system.time_step *= 1.
        if espressomd.has_features("ELECTROSTATICS"):
            self.system.electrostatics.solver = espressomd.electrostatics.DH(
                prefactor=1., kappa=1., r_cut=1.)  # should not fail
            self.system.electrostatics.clear()
        with self.assertRaisesRegex(RuntimeError, "MD cell geometry change not supported by EK"):
            self.system.box_l = [1., 2., 3.]
        np.testing.assert_allclose(np.copy(self.system.box_l), 6., atol=1e-7)

    def test_ek_reactants(self):
        ek_species = self.make_default_ek_species()
        ek_reactant = espressomd.electrokinetics.EKReactant(
            ekspecies=ek_species, stoech_coeff=-2.0, order=2.0)
        self.assertAlmostEqual(ek_reactant.stoech_coeff, -2.0, delta=self.atol)
        self.assertAlmostEqual(ek_reactant.order, 2.0, delta=self.atol)
        ek_reactant.stoech_coeff = 1.0
        self.assertAlmostEqual(ek_reactant.stoech_coeff, 1.0, delta=self.atol)

        with self.assertRaisesRegex(RuntimeError, f"(Parameter|Property) 'order' is read-only"):
            ek_reactant.order = 1.5

    def test_ek_indexed_reactions(self):
        ek_species = self.make_default_ek_species()
        ek_reactant = espressomd.electrokinetics.EKReactant(
            ekspecies=ek_species, stoech_coeff=-2.0, order=2.0)
        ek_reaction = espressomd.electrokinetics.EKIndexedReaction(
            reactants=[ek_reactant], coefficient=1.5, lattice=self.lattice, tau=self.params["tau"])
        self.assertAlmostEqual(ek_reaction.coefficient, 1.5, delta=self.atol)
        ek_reaction.coefficient = 0.5
        self.assertAlmostEqual(ek_reaction.coefficient, 0.5, delta=self.atol)

        # -- node interface --
        self.assertFalse(ek_reaction[1, 1, 1])
        ek_reaction[1, 1, 1] = True
        self.assertTrue(ek_reaction[1, 1, 1])
        ek_reaction.remove_node_from_index([1, 1, 1])
        self.assertFalse(ek_reaction[1, 1, 1])
        ek_reaction.add_node_to_index([1, 1, 1])
        self.assertTrue(ek_reaction[1, 1, 1])

        # multiple independent nodes can be set/cleared independently
        ek_reaction[2, 3, 4] = True
        self.assertTrue(ek_reaction[1, 1, 1])
        self.assertTrue(ek_reaction[2, 3, 4])
        self.assertFalse(ek_reaction[0, 0, 0])
        ek_reaction[1, 1, 1] = False
        self.assertFalse(ek_reaction[1, 1, 1])
        self.assertTrue(ek_reaction[2, 3, 4])

        # bad key type and bad key length both raise TypeError
        with self.assertRaisesRegex(TypeError, "is not a valid index"):
            ek_reaction["bad"]
        with self.assertRaisesRegex(TypeError, "is not a valid index"):
            ek_reaction[0, 1]

        # -- slice interface --
        nx, ny, nz = ek_reaction.shape

        # initial state: all nodes are out of the reaction index
        is_boundary = ek_reaction[:, :, :].is_boundary
        self.assertEqual(is_boundary.shape, (nx, ny, nz))
        self.assertEqual(is_boundary.dtype, bool)
        # a node set via the node interface shows up in a slice read
        np.testing.assert_array_equal(is_boundary[2, 3, 4], True)
        # clear it and verify the full grid is False
        ek_reaction[:, :, :] = False
        np.testing.assert_array_equal(ek_reaction[:, :, :].is_boundary, False)

        # returned array is locked (non-writable)
        locked = ek_reaction[:, :, :].is_boundary
        with self.assertRaisesRegex(ValueError,
                                    "ESPResSo array properties return non-writable arrays"):
            locked[0, 0, 0] = True

        # set a 2D slice (integer dim collapses from the shape) and read back
        ek_reaction[1, :, :] = True
        result = ek_reaction[1, :, :].is_boundary
        self.assertEqual(result.shape, (ny, nz))
        np.testing.assert_array_equal(result, True)
        # cells outside the written slice remain False
        np.testing.assert_array_equal(ek_reaction[0, :, :].is_boundary, False)
        np.testing.assert_array_equal(ek_reaction[2:, :, :].is_boundary, False)

        # reset, then broadcast a scalar to a 3D slice
        ek_reaction[:, :, :] = False
        ek_reaction[3:5, :, :] = True
        result = ek_reaction[3:5, :, :].is_boundary
        self.assertEqual(result.shape, (2, ny, nz))
        np.testing.assert_array_equal(result, True)

        # write a shaped bool array and read it back (round-trip)
        ek_reaction[:, :, :] = False
        values = np.zeros((2, ny, nz), dtype=bool)
        values[0, 0, 1] = True
        values[1, ny - 1, 0] = True
        ek_reaction[3:5, :, :] = values
        np.testing.assert_array_equal(
            ek_reaction[3:5, :, :].is_boundary, values)

        # getters
        ek_slice = ek_reaction[0:2, -4:-1, 1:2]
        np.testing.assert_array_equal(
            np.copy(ek_slice.call_method("get_slice_size")), [2, 3, 1])
        np.testing.assert_array_equal(
            np.copy(ek_slice.call_method("get_slice_ranges")),
            [[0, 8, 1], [2, 11, 2]])
        self.assertEqual(ek_slice.call_method("get_reaction_sip"), ek_reaction)

        # wrong shape raises ValueError
        with self.assertRaisesRegex(ValueError,
                                    "Input-dimensions of 'is_boundary' array"):
            ek_reaction[1, :, :] = np.zeros((ny + 1, nz), dtype=bool)

        # empty (out-of-bounds) slice: getter returns an empty array
        empty_slice = ek_reaction[nx:nx + 5, :, :]
        self.assertEqual(empty_slice.is_boundary.shape[0], 0)

        # empty slice setter raises AttributeError
        with self.assertRaisesRegex(
                AttributeError,
                "Cannot set properties of an empty 'EKIndexedReactionSlice' object"):
            empty_slice.is_boundary = True

        # unknown property name raises RuntimeError
        with self.assertRaisesRegex(RuntimeError,
                                    "Unknown EK indexed reaction property 'unknown'"):
            ek_reaction[:, :, :].call_method("get_value_shape", name="unknown")

    def test_ek_fluctuations(self):
        """
        smoke test, see `ek_fluctuations.py` for a statistical test
        """
        lattice = espressomd.electrokinetics.Lattice(
            n_ghost_layers=1, agrid=1.)
        ek_solver = espressomd.electrokinetics.EKNone(
            lattice=lattice, tau=self.params["tau"], **self.ek_params)
        ek_species = self.ek_species_class(
            lattice=lattice, density=20., valency=0., advection=False,
            diffusion=0.01, friction_coupling=False, thermalized=True, seed=42,
            tau=self.params["tau"], **self.ek_params)
        self.assertTrue(ek_species.thermalized)
        self.assertEqual(ek_species.seed, 42)
        self.system.ekcontainer.solver = ek_solver
        self.system.ekcontainer.add(ek_species)
        self.system.integrator.run(20)

    @utx.skipIfMissingFeatures(["WALBERLA_FFT"])
    def test_ek_friction_advection(self):
        """
        smoke test, see `ek_eof.py` for a statistical test
        """
        lattice = espressomd.electrokinetics.Lattice(
            n_ghost_layers=1, agrid=1.)
        lb_fluid = self.lb_fluid_class(
            lattice=lattice, density=1., kinematic_viscosity=1. / 6.,
            tau=self.params["tau"], **self.lb_params)
        ek_species = self.ek_species_class(
            lattice=lattice, density=1., kT=2., valency=1.1, diffusion=0.25,
            friction_coupling=True, advection=True, ext_efield=[0., 0.001, 0.],
            tau=self.params["tau"], **self.ek_params)
        ek_wallcharge = self.ek_species_class(
            lattice=lattice, density=0., kT=2., valency=-1.1, diffusion=0.,
            friction_coupling=False, advection=False, ext_efield=[0., 0., 0.],
            tau=self.params["tau"], **self.ek_params)
        ek_solver = self.ek_solver_class(
            lattice=lattice, permittivity=0.3, tau=self.params["tau"],
            **self.ek_params)
        self.assertTrue(ek_species.friction_coupling)
        self.assertTrue(ek_species.advection)
        self.assertFalse(ek_wallcharge.friction_coupling)
        self.assertFalse(ek_wallcharge.advection)
        self.assertAlmostEqual(ek_solver.permittivity, 0.3, delta=self.atol)
        self.system.ekcontainer.solver = ek_solver
        self.system.ekcontainer.add(ek_species)
        self.system.ekcontainer.add(ek_wallcharge)
        self.system.lb = lb_fluid
        self.system.integrator.run(20)

    def test_grid_index(self):
        ek_species = self.make_default_ek_species()
        ek_solver = self.system.ekcontainer.solver
        ek_reactant = espressomd.electrokinetics.EKReactant(
            ekspecies=ek_species, stoech_coeff=-2.0, order=2.0)
        ek_reaction = espressomd.electrokinetics.EKIndexedReaction(
            reactants=[ek_reactant], coefficient=1.5, lattice=self.lattice, tau=self.params["tau"])
        # check ranges and out-of-bounds access
        shape = np.around(self.system.box_l / self.params["agrid"]).astype(int)
        int_shape = [int(x) for x in shape]  # cast away numpy integer types
        for i in range(3):
            n = [0, 0, 0]
            n[i] -= shape[i]
            ek_reaction[n[0], n[1], n[2]] = True
            self.assertTrue(ek_reaction[0, 0, 0])
            self.assertEqual(ek_reaction[tuple(n)], ek_reaction[0, 0, 0])
            self.assertEqual(ek_species[tuple(n)], ek_species[0, 0, 0])
            self.assertEqual(ek_solver[tuple(n)], ek_solver[0, 0, 0])
            for offset in (int_shape[i] + 1, -(int_shape[i] + 1)):
                n = [0, 0, 0]
                n[i] += offset
                err_msg = rf"provided index \[{str(n)[1:-1]}\] is out of range for shape \[{str(int_shape)[1:-1]}\]"  # nopep8
                with self.assertRaisesRegex(IndexError, err_msg):
                    ek_reaction[tuple(n)]
                with self.assertRaisesRegex(IndexError, err_msg):
                    ek_species[tuple(n)]
                with self.assertRaisesRegex(IndexError, err_msg):
                    ek_solver[tuple(n)]
        # node index
        for ek_obj in [ek_species, ek_solver]:
            node = ek_obj[1, 2, 3]
            with self.assertRaisesRegex(RuntimeError, "Parameter 'index' is read-only"):
                node.index = [2, 4, 6]
            np.testing.assert_array_equal(node.index, [1, 2, 3])
            retval = node.call_method("override_index", index=[2, 4, 6])
            self.assertEqual(retval, 0)
            np.testing.assert_array_equal(node.index, [2, 4, 6])
            retval = node.call_method("override_index", index=[0, 0, shape[2]])
            self.assertEqual(retval, 1)
            np.testing.assert_array_equal(node.index, [2, 4, 6])
            np.testing.assert_array_equal(ek_obj[-1, -1, -1].index, shape - 1)

    def test_runtime_exceptions(self):
        # set up a valid species
        ek_species = self.make_default_ek_species()
        ek_species.kT = 0.
        self.system.ekcontainer.add(ek_species)
        self.system.integrator.run(1)

        # check exceptions without LB force field
        with self.assertRaisesRegex(Exception, "friction coupling enabled but no force field accessible"):
            ek_species.friction_coupling = True
            ek_species.advection = False
            self.system.integrator.run(1)

        # check exceptions without LB velocity field
        with self.assertRaisesRegex(Exception, "advection enabled but no velocity field accessible"):
            ek_species.friction_coupling = False
            ek_species.advection = True
            self.system.integrator.run(1)

        # non-diffusive species don't trigger exceptions due to early exit
        ek_species.friction_coupling = True
        ek_species.advection = True
        ek_species.diffusion = 0.
        self.system.integrator.run(1)

        # check exceptions with an incompatible LB time step
        with self.assertRaisesRegex(Exception, "LB and EK are active but with different time steps"):
            lb = self.lb_fluid_class(
                lattice=self.lattice, density=0.5, kinematic_viscosity=3.,
                tau=2. * self.params["tau"], **self.lb_params)
            self.system.lb = lb
            self.system.integrator.run(1)

        # reset global variable fluid_step
        self.system.ekcontainer.clear()
        self.system.integrator.run(1)

    def test_ek_bulk_reactions(self):
        ek_species = self.make_default_ek_species()
        ek_reactant = espressomd.electrokinetics.EKReactant(
            ekspecies=ek_species, stoech_coeff=-2.0, order=2.0)
        ek_reaction = espressomd.electrokinetics.EKBulkReaction(
            reactants=[ek_reactant], coefficient=1.5, lattice=self.lattice, tau=self.params["tau"])
        self.assertAlmostEqual(ek_reaction.coefficient, 1.5, delta=self.atol)
        ek_reaction.coefficient = 0.5
        self.assertAlmostEqual(ek_reaction.coefficient, 0.5, delta=self.atol)

    def test_raise_if_read_only(self):
        ek_species = self.make_default_ek_species()
        for key in {"lattice", "shape",
                    "single_precision", "seed", "thermalized"}:
            with self.assertRaisesRegex(RuntimeError, f"(Parameter|Property) '{key}' is read-only"):
                setattr(ek_species, key, 0)

    def test_ctor_exceptions(self):
        def make_kwargs(**kwargs):
            ek_kwargs = {"lattice": self.lattice}
            ek_kwargs.update(self.ek_species_params)
            ek_kwargs.update(self.ek_params)
            ek_kwargs.update(kwargs)
            return ek_kwargs

        with self.assertRaisesRegex(ValueError, "Parameter 'tau' must be > 0"):
            self.ek_species_class(**make_kwargs(tau=0.))
        with self.assertRaisesRegex(ValueError, "Parameter 'density' must be >= 0"):
            self.ek_species_class(**make_kwargs(density=-1.))
        with self.assertRaisesRegex(ValueError, "Parameter 'kT' must be >= 0"):
            self.ek_species_class(**make_kwargs(kT=-1.))
        with self.assertRaisesRegex(ValueError, "Parameter 'tau' must be > 0"):
            espressomd.electrokinetics.EKContainer(
                tau=0., solver=self.system.ekcontainer.solver)
        for thermalized in (True, False):
            with self.assertRaisesRegex(ValueError, "Parameter 'seed' must be >= 0"):
                self.ek_species_class(
                    **make_kwargs(thermalized=thermalized, seed=-1))
        with self.assertRaisesRegex(ValueError, "Parameter 'seed' is required for thermalized EKSpecies"):
            self.ek_species_class(**make_kwargs(thermalized=True))

        # when ekcontainer is None, no solver can be attached
        self.system.ekcontainer = None
        self.system.ekcontainer.solver = None
        with self.assertRaisesRegex(RuntimeError, "Parameter 'solver' is read-only"):
            self.system.ekcontainer.solver = espressomd.electrokinetics.EKNone(
                lattice=self.lattice, tau=self.params["tau"], **self.ek_params)
        self.assertIsNone(self.system.ekcontainer.solver)

    def test_rollback(self):
        """check rollback to a valid state when setter fails"""
        node_grid = np.copy(self.system.cell_system.node_grid)
        world_size = np.prod(node_grid)
        if world_size <= 4:
            wrong_box_l = [1., 1., 7.] if world_size == 1 else 2. * node_grid
            lattice1 = espressomd.electrokinetics.Lattice(
                n_ghost_layers=2, agrid=1., box_l=self.system.box_l)
            lattice2 = espressomd.electrokinetics.Lattice(
                n_ghost_layers=2, agrid=1., box_l=wrong_box_l)
            solver_valid = espressomd.electrokinetics.EKNone(
                lattice=lattice1, tau=self.params["tau"], **self.ek_params)
            solver_wrong = espressomd.electrokinetics.EKNone(
                lattice=lattice2, tau=self.params["tau"], **self.ek_params)
            self.system.ekcontainer = espressomd.electrokinetics.EKContainer(
                tau=self.system.time_step, solver=solver_valid)
            with self.assertRaisesRegex(RuntimeError, "waLBerla and ESPResSo disagree about domain decomposition"):
                self.system.ekcontainer = espressomd.electrokinetics.EKContainer(
                    tau=self.system.time_step, solver=solver_wrong)
            self.assertEqual(self.system.ekcontainer.solver, solver_valid)

    def test_node_grid_change(self):
        """check MPI Cartesian communicator invalidation"""
        node_grid = np.copy(self.system.cell_system.node_grid)
        # create a species, slice and node for the current MPI topology
        ek_solver = espressomd.electrokinetics.EKNone(
            lattice=self.lattice, tau=self.params["tau"], **self.ek_params)
        ek_species = self.make_default_ek_species()
        ek_node = ek_species[0, 0, 0]
        ek_slice = ek_species[0:5, 0:5, 0:5]
        self.system.ekcontainer = espressomd.electrokinetics.EKContainer(
            tau=self.system.time_step, solver=ek_solver)
        self.system.ekcontainer.add(ek_species)
        # veto node grid change
        with self.assertRaisesRegex(RuntimeError, "MPI topology change not supported by EK"):
            self.system.cell_system.node_grid = node_grid
        self.system.ekcontainer = None
        # invalidate MPI Cartesian communicator
        self.system.cell_system.node_grid = node_grid
        # create a new species
        ek_solver_new = espressomd.electrokinetics.EKNone(
            lattice=self.lattice, tau=self.params["tau"], **self.ek_params)
        ek_species_new = self.make_default_ek_species()
        self.system.ekcontainer = espressomd.electrokinetics.EKContainer(
            tau=self.system.time_step, solver=ek_solver_new)
        self.system.ekcontainer.add(ek_species_new)
        # prevent binding of an expired EK object
        with self.assertRaisesRegex(RuntimeError, "the MPI Cartesian communicator of this EK object has expired"):
            self.system.ekcontainer.add(ek_species)
        self.assertEqual(len(self.system.ekcontainer), 1)
        self.assertEqual(self.system.ekcontainer[0], ek_species_new)
        # expired MPI communicator doesn't prevent read access to the fields
        _ = ek_node.density
        _ = ek_slice.density
        # expired MPI communicator prevents write access to the fields
        for handle in [ek_node, ek_slice,
                       ek_species[0, 0, 0], ek_species[0:5, 0:5, 0:5]]:
            with self.assertRaisesRegex(RuntimeError, "the MPI Cartesian communicator of this EK object has expired"):
                handle.density = 1.

    def test_bool_operations_on_node(self):
        ekspecies = self.make_default_ek_species()
        # test __eq()__ where a node is equal to itself and not equal to any
        # other node
        assert ekspecies[0, 0, 0] == ekspecies[0, 0, 0]
        shape = np.around(self.system.box_l / self.params["agrid"]).astype(int)
        nodes = [
            ekspecies[ijk] for ijk in itertools.product(
                range(shape[0]), range(shape[1]), range(shape[2]))]
        nodes.remove(ekspecies[0, 0, 0])
        assert all(ekspecies[0, 0, 0] != node for node in nodes)
        # test __hash()__ intercept to identify nodes based on index rather
        # than name. set() constructor runs hash()
        subset1, subset2 = nodes[:-10], nodes[-10:]
        assert len(set(subset1 + subset1)) == len(subset1)
        assert len(set(subset1 + subset2)) == len(subset1) + len(subset2)


@utx.skipIfMissingFeatures(["WALBERLA"])
class EKTestWalberlaDoublePrecisionCPU(EKTest, ut.TestCase):

    """Test for the waLBerla implementation of the EK in double-precision."""

    lb_fluid_class = espressomd.lb.LBFluid
    ek_lattice_class = espressomd.electrokinetics.Lattice
    ek_species_class = espressomd.electrokinetics.EKSpecies
    ek_solver_class = espressomd.electrokinetics.EKFFT
    ek_params = {"single_precision": False, "gpu": False}
    lb_params = {"single_precision": False, "gpu": False}
    atol = 1e-10
    rtol = 1e-7


@utx.skipIfMissingFeatures(["WALBERLA"])
class EKTestWalberlaSinglePrecisionCPU(EKTest, ut.TestCase):

    """Test for the waLBerla implementation of the EK in single-precision."""

    lb_fluid_class = espressomd.lb.LBFluid
    ek_lattice_class = espressomd.electrokinetics.Lattice
    ek_species_class = espressomd.electrokinetics.EKSpecies
    ek_solver_class = espressomd.electrokinetics.EKFFT
    ek_params = {"single_precision": True, "gpu": False}
    lb_params = {"single_precision": True, "gpu": False}
    atol = 1e-7
    rtol = 5e-5


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class EKTestWalberlaDoublePrecisionGPU(EKTest, ut.TestCase):

    """Test for the waLBerla implementation of the EK in double-precision."""

    lb_fluid_class = espressomd.lb.LBFluid
    ek_lattice_class = espressomd.electrokinetics.Lattice
    ek_species_class = espressomd.electrokinetics.EKSpecies
    ek_solver_class = espressomd.electrokinetics.EKFFT
    ek_params = {"single_precision": False, "gpu": True}
    lb_params = {"single_precision": False, "gpu": True}
    atol = 1e-10
    rtol = 1e-7


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class EKTestWalberlaSinglePrecisionGPU(EKTest, ut.TestCase):

    """Test for the waLBerla implementation of the EK in single-precision."""

    lb_fluid_class = espressomd.lb.LBFluid
    ek_lattice_class = espressomd.electrokinetics.Lattice
    ek_species_class = espressomd.electrokinetics.EKSpecies
    ek_solver_class = espressomd.electrokinetics.EKFFT
    ek_params = {"single_precision": True, "gpu": True}
    lb_params = {"single_precision": True, "gpu": True}
    atol = 1e-7
    rtol = 5e-5


if __name__ == "__main__":
    ut.main()
