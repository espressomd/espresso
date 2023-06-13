#
# Copyright (C) 2010-2023 The ESPResSo project
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

import sys
import numpy as np
import itertools
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
        self.lattice = self.ek_lattice_class(
            n_ghost_layers=1, agrid=self.params["agrid"])
        ek_solver = espressomd.electrokinetics.EKNone(lattice=self.lattice)
        self.system.ekcontainer.solver = ek_solver
        self.system.ekcontainer.tau = self.system.time_step

    def tearDown(self):
        self.system.actors.clear()
        self.system.ekcontainer.clear()
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.time_step = self.params["tau"]
        self.system.ekcontainer.solver = None

    def make_default_ek_species(self):
        return self.ek_species_class(
            lattice=self.lattice, **self.ek_params, **self.ek_species_params)

    def test_ek_container_poisson_solver(self):
        ekcontainer = self.system.ekcontainer
        ekcontainer.solver = None
        self.assertFalse(ekcontainer.call_method("is_poisson_solver_set"))
        self.assertIsNone(ekcontainer.solver)
        ek_solver = espressomd.electrokinetics.EKNone(lattice=self.lattice)
        self.assertFalse(ekcontainer.call_method("is_poisson_solver_set"))
        ekcontainer.solver = ek_solver
        self.assertTrue(ekcontainer.call_method("is_poisson_solver_set"))
        self.assertIsInstance(self.system.ekcontainer.solver,
                              espressomd.electrokinetics.EKNone)

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

    def check_ek_species_properties(self, species):
        agrid = self.params["agrid"]
        # check getters
        self.assertEqual(species.lattice.n_ghost_layers, 1)
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
        ek_solver = espressomd.electrokinetics.EKFFT(
            lattice=self.lattice, permittivity=0.01,
            single_precision=self.ek_params["single_precision"])
        self.assertEqual(
            ek_solver.single_precision,
            self.ek_params["single_precision"])
        self.assertAlmostEqual(ek_solver.permittivity, 0.01, delta=self.atol)
        ek_solver.permittivity = 0.05
        self.assertAlmostEqual(ek_solver.permittivity, 0.05, delta=self.atol)

        self.system.ekcontainer.solver = ek_solver
        self.assertTrue(
            self.system.ekcontainer.call_method("is_poisson_solver_set"))
        self.assertIsInstance(self.system.ekcontainer.solver,
                              espressomd.electrokinetics.EKFFT)

    def test_ek_none_solver(self):
        ek_solver = espressomd.electrokinetics.EKNone(
            lattice=self.lattice,
            single_precision=self.ek_params["single_precision"])
        self.assertEqual(
            ek_solver.single_precision,
            self.ek_params["single_precision"])

        self.system.ekcontainer.solver = ek_solver
        self.assertTrue(
            self.system.ekcontainer.call_method("is_poisson_solver_set"))
        self.assertIsInstance(self.system.ekcontainer.solver,
                              espressomd.electrokinetics.EKNone)

    def test_ek_solver_exceptions(self):
        ek_solver = self.system.ekcontainer.solver
        ek_species = self.make_default_ek_species()
        self.system.ekcontainer.add(ek_species)
        incompatible_lattice = self.ek_lattice_class(
            n_ghost_layers=2, agrid=self.params["agrid"])
        incompatible_ek_species = self.ek_species_class(
            lattice=incompatible_lattice, **self.ek_params, **self.ek_species_params)
        with self.assertRaisesRegex(RuntimeError, "EKSpecies lattice incompatible with existing Poisson solver lattice"):
            self.system.ekcontainer.add(incompatible_ek_species)
        with self.assertRaisesRegex(RuntimeError, "EKSpecies lattice incompatible with existing EKSpecies lattice"):
            self.system.ekcontainer.solver = None
            self.system.ekcontainer.add(incompatible_ek_species)
        with self.assertRaisesRegex(ValueError, "Parameter 'tau' is required when container isn't empty"):
            self.system.ekcontainer.tau = None
        with self.assertRaisesRegex(RuntimeError, "Poisson solver lattice incompatible with existing EKSpecies lattice"):
            self.system.ekcontainer.clear()
            self.system.ekcontainer.add(incompatible_ek_species)
            self.system.ekcontainer.solver = ek_solver
        with self.assertRaisesRegex(ValueError, "Parameter 'tau' must be > 0"):
            self.system.ekcontainer.tau = 0.
        self.assertAlmostEqual(
            self.system.ekcontainer.tau, self.system.time_step, delta=1e-7)
        self.system.ekcontainer.clear()
        self.system.ekcontainer.tau = None
        self.assertIsNone(self.system.ekcontainer.tau)

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
        # boundaries
        self.assertFalse(ek_reaction[1, 1, 1])
        ek_reaction[1, 1, 1] = True
        self.assertTrue(ek_reaction[1, 1, 1])
        ek_reaction.remove_node_from_index([1, 1, 1])
        self.assertFalse(ek_reaction[1, 1, 1])
        ek_reaction.add_node_to_index([1, 1, 1])
        self.assertTrue(ek_reaction[1, 1, 1])

    def test_grid_index(self):
        ek_species = self.make_default_ek_species()
        ek_reactant = espressomd.electrokinetics.EKReactant(
            ekspecies=ek_species, stoech_coeff=-2.0, order=2.0)
        ek_reaction = espressomd.electrokinetics.EKIndexedReaction(
            reactants=[ek_reactant], coefficient=1.5, lattice=self.lattice, tau=self.params["tau"])
        # check ranges and out-of-bounds access
        shape = np.around(self.system.box_l / self.params["agrid"]).astype(int)
        for i in range(3):
            n = [0, 0, 0]
            n[i] -= shape[i]
            ek_reaction[n[0], n[1], n[2]] = True
            self.assertTrue(ek_reaction[0, 0, 0])
            self.assertEqual(ek_reaction[tuple(n)], ek_reaction[0, 0, 0])
            self.assertEqual(ek_species[tuple(n)], ek_species[0, 0, 0])
            for offset in (shape[i] + 1, -(shape[i] + 1)):
                n = [0, 0, 0]
                n[i] += offset
                err_msg = rf"provided index \[{str(n)[1:-1]}\] is out of range for shape \[{str(list(shape))[1:-1]}\]"
                with self.assertRaisesRegex(IndexError, err_msg):
                    ek_reaction[tuple(n)]
                with self.assertRaisesRegex(IndexError, err_msg):
                    ek_species[tuple(n)]
        # node index
        node = ek_species[1, 2, 3]
        with self.assertRaisesRegex(RuntimeError, "Parameter 'index' is read-only"):
            node.index = [2, 4, 6]
        np.testing.assert_array_equal(node.index, [1, 2, 3])
        retval = node.call_method("override_index", index=[2, 4, 6])
        self.assertEqual(retval, 0)
        np.testing.assert_array_equal(node.index, [2, 4, 6])
        retval = node.call_method("override_index", index=[0, 0, shape[2]])
        self.assertEqual(retval, 1)
        np.testing.assert_array_equal(node.index, [2, 4, 6])
        np.testing.assert_array_equal(ek_species[-1, -1, -1].index, shape - 1)

    def test_runtime_exceptions(self):
        # set up a valid species
        ek_species = self.make_default_ek_species()
        ek_species.kT = 0.
        self.system.ekcontainer.add(ek_species)
        self.system.integrator.run(1)

        print("\nTesting EK runtime error messages:", file=sys.stderr)
        sys.stderr.flush()

        # check exceptions without EK Poisson solver
        with self.assertRaisesRegex(Exception, "EK requires a Poisson solver"):
            self.system.ekcontainer.solver = None
            self.system.integrator.run(1)
        ek_solver = espressomd.electrokinetics.EKNone(lattice=self.lattice)
        self.system.ekcontainer.solver = ek_solver

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
            self.system.actors.add(lb)
            self.system.integrator.run(1)

        print("End of EK runtime error messages", file=sys.stderr)
        sys.stderr.flush()

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
        for key in {"lattice", "shape", "single_precision"}:
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
class EKTestWalberla(EKTest, ut.TestCase):

    """Test for the Walberla implementation of the EK in double-precision."""

    lb_fluid_class = espressomd.lb.LBFluidWalberla
    ek_lattice_class = espressomd.electrokinetics.LatticeWalberla
    ek_species_class = espressomd.electrokinetics.EKSpecies
    ek_params = {"single_precision": False}
    lb_params = {"single_precision": False}
    atol = 1e-10
    rtol = 1e-7


@utx.skipIfMissingFeatures(["WALBERLA"])
class EKTestWalberlaSinglePrecision(EKTest, ut.TestCase):

    """Test for the Walberla implementation of the EK in single-precision."""

    lb_fluid_class = espressomd.lb.LBFluidWalberla
    ek_lattice_class = espressomd.electrokinetics.LatticeWalberla
    ek_species_class = espressomd.electrokinetics.EKSpecies
    ek_params = {"single_precision": True}
    lb_params = {"single_precision": True}
    atol = 1e-7
    rtol = 5e-5


if __name__ == "__main__":
    ut.main()
