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
import numpy as np
import itertools
import sys

import espressomd
import espressomd.lb
import espressomd.EKSpecies


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
                         "ext_efield": [0.1, 0.2, 0.3]}

    system.periodicity = [True, True, True]
    system.time_step = params["tau"]
    system.cell_system.skin = 1.0
    ek_solver_set = False

    @classmethod
    def setUpClass(cls):
        cls.lattice = cls.ek_lattice_class(
            n_ghost_layers=1, agrid=cls.params["agrid"])

    def tearDown(self):
        self.system.actors.clear()
        self.system.ekcontainer.clear()
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.time_step = self.params["tau"]

    def make_default_ek_species(self):
        return self.ek_species_class(
            lattice=self.lattice,
            single_precision=self.ek_params["single_precision"],
            **self.ek_species_params)

    def test_00_ek_container(self):
        def check():
            status = self.system.ekcontainer.call_method(
                "is_poissonsolver_set")
            self.assertEqual(status, EKTest.ek_solver_set)

        check()
        ek_solver = espressomd.EKSpecies.EKNone(lattice=self.lattice)
        check()
        self.system.ekcontainer.solver = ek_solver
        EKTest.ek_solver_set = True
        check()

    def test_ek_species(self):
        # inactive species
        ek_species = self.make_default_ek_species()
        self.check_ek_species_properties(ek_species)

        ek_solver = espressomd.EKSpecies.EKNone(lattice=self.lattice)
        self.system.ekcontainer.tau = self.system.time_step
        self.system.ekcontainer.solver = ek_solver
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
        espressomd.EKSpecies.EKReactant(
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
            species.is_single_precision,
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
        node.flux_boundary = espressomd.EKSpecies.FluxBoundary([1., 2., 3.])
        self.assertIsInstance(
            node.flux_boundary,
            espressomd.EKSpecies.FluxBoundary)
        np.testing.assert_allclose(
            np.copy(node.flux_boundary.flux), [1., 2., 3.], atol=self.atol)
        node.density_boundary = espressomd.EKSpecies.DensityBoundary(4.5)
        self.assertIsInstance(
            node.density_boundary,
            espressomd.EKSpecies.DensityBoundary)
        np.testing.assert_allclose(
            np.copy(node.density_boundary.density), 4.5, atol=self.atol)
        node.density_boundary = None
        node.flux_boundary = None
        self.assertIsNone(node.density_boundary)
        self.assertIsNone(node.flux_boundary)
        with self.assertRaisesRegex(TypeError, "must be an instance of DensityBoundary or None"):
            node.density_boundary = 4.6
        with self.assertRaisesRegex(TypeError, "must be an instance of FluxBoundary or None"):
            node.flux_boundary = 4.6

    def test_ek_fft_solvers(self):
        ek_solver = espressomd.EKSpecies.EKFFT(
            lattice=self.lattice, permittivity=0.01,
            single_precision=self.ek_params["single_precision"])
        self.assertEqual(
            ek_solver.is_single_precision,
            self.ek_params["single_precision"])
        self.assertAlmostEqual(ek_solver.permittivity, 0.01, delta=self.atol)
        ek_solver.permittivity = 0.05
        self.assertAlmostEqual(ek_solver.permittivity, 0.05, delta=self.atol)

    def test_ek_none_solvers(self):
        ek_solver = espressomd.EKSpecies.EKNone(
            lattice=self.lattice,
            single_precision=self.ek_params["single_precision"])
        self.assertEqual(
            ek_solver.is_single_precision,
            self.ek_params["single_precision"])

    def test_ek_reactants(self):
        ek_species = self.make_default_ek_species()
        ek_reactant = espressomd.EKSpecies.EKReactant(
            ekspecies=ek_species, stoech_coeff=-2.0, order=2.0)
        self.assertAlmostEqual(ek_reactant.stoech_coeff, -2.0, delta=self.atol)
        self.assertAlmostEqual(ek_reactant.order, 2.0, delta=self.atol)
        ek_reactant.stoech_coeff = 1.0
        ek_reactant.order = 1.5
        self.assertAlmostEqual(ek_reactant.stoech_coeff, 1.0, delta=self.atol)
        self.assertAlmostEqual(ek_reactant.order, 1.5, delta=self.atol)

    def test_ek_indexed_reactions(self):
        ek_species = self.make_default_ek_species()
        ek_reactant = espressomd.EKSpecies.EKReactant(
            ekspecies=ek_species, stoech_coeff=-2.0, order=2.0)
        ek_reaction = espressomd.EKSpecies.EKIndexedReaction(
            reactants=[ek_reactant], coefficient=1.5, lattice=self.lattice)
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
        ek_reactant = espressomd.EKSpecies.EKReactant(
            ekspecies=ek_species, stoech_coeff=-2.0, order=2.0)
        ek_reaction = espressomd.EKSpecies.EKIndexedReaction(
            reactants=[ek_reactant], coefficient=1.5, lattice=self.lattice)
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
        with self.assertRaisesRegex(Exception, "LB and EK are active but with different timesteps"):
            lb = self.lb_fluid_class(
                lattice=self.lattice, density=0.5, viscosity=3.,
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
        ek_reactant = espressomd.EKSpecies.EKReactant(
            ekspecies=ek_species, stoech_coeff=-2.0, order=2.0)
        ek_reaction = espressomd.EKSpecies.EKBulkReaction(
            reactants=[ek_reactant], coefficient=1.5, lattice=self.lattice)
        self.assertAlmostEqual(ek_reaction.coefficient, 1.5, delta=self.atol)
        ek_reaction.coefficient = 0.5
        self.assertAlmostEqual(ek_reaction.coefficient, 0.5, delta=self.atol)

    def test_raise_if_read_only(self):
        ek_species = self.make_default_ek_species()
        for key in {"lattice", "shape", "is_single_precision"}:
            with self.assertRaisesRegex(RuntimeError, f"(Parameter|Property) '{key}' is read-only"):
                setattr(ek_species, key, 0)

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


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class EKTestWalberla(EKTest, ut.TestCase):

    """Test for the Walberla implementation of the EK in double-precision."""

    lb_fluid_class = espressomd.lb.LBFluidWalberla
    ek_lattice_class = espressomd.lb.LatticeWalberla
    ek_species_class = espressomd.EKSpecies.EKSpecies
    ek_params = {"single_precision": False}
    lb_params = {"single_precision": False}
    atol = 1e-10
    rtol = 1e-7


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class EKTestWalberlaSinglePrecision(EKTest, ut.TestCase):

    """Test for the Walberla implementation of the EK in single-precision."""

    lb_fluid_class = espressomd.lb.LBFluidWalberla
    ek_lattice_class = espressomd.lb.LatticeWalberla
    ek_species_class = espressomd.EKSpecies.EKSpecies
    ek_params = {"single_precision": True}
    lb_params = {"single_precision": True}
    atol = 1e-7
    rtol = 5e-5


if __name__ == "__main__":
    ut.main()
