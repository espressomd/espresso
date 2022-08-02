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
                         "ext_efield": [0., 0., 0.]}

    system.periodicity = [True, True, True]
    system.time_step = params["tau"]
    system.cell_system.skin = 1.0

    @classmethod
    def setUpClass(cls):
        cls.lattice = cls.ek_lattice_class(
            n_ghost_layers=1, agrid=cls.params["agrid"])

    def tearDown(self):
        self.system.ekcontainer.clear()
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.time_step = self.params["tau"]

    def make_default_ek_species(self):
        return self.ek_species_class(
            lattice=self.lattice,
            single_precision=self.ek_params["single_precision"],
            **self.ek_species_params)

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
        # check EK object
        self.assertEqual(species.lattice.n_ghost_layers, 1)
        self.assertAlmostEqual(species.lattice.agrid, agrid, delta=self.atol)
        # self.assertAlmostEqual(species.density, 0.85, delta=self.atol) # TODO
        self.assertAlmostEqual(species.diffusion, 0.1, delta=self.atol)
        self.assertAlmostEqual(species.valency, 0.0, delta=self.atol)
        self.assertAlmostEqual(species.kT, 1.5, delta=self.atol)
        self.assertFalse(species.advection)
        self.assertFalse(species.friction_coupling)
        self.assertEqual(
            species.is_single_precision,
            self.ek_params["single_precision"])
        np.testing.assert_allclose(
            np.copy(species.ext_efield), [0., 0., 0.], atol=self.atol)
        species.diffusion = 0.2
        self.assertAlmostEqual(species.diffusion, 0.2, delta=self.atol)
        ext_f = [0.01, 0.02, 0.03]
        species.ext_efield = ext_f
        np.testing.assert_allclose(
            np.copy(species.ext_efield), ext_f, atol=self.atol)
        # check node getters/setters
        self.assertAlmostEqual(species[0, 0, 0].density, 0.85, delta=self.atol)
        species[0, 0, 0].density = 0.90
        self.assertAlmostEqual(species[0, 0, 0].density, 0.90, delta=self.atol)
        with self.assertRaises(RuntimeError):
            species[0, 0, 0].density = [1, 2]
        with self.assertRaises(TypeError):
            species[0, 1].density = 1.

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

    # TODO walberla: fix infinite loop
#    def test_grid_index(self):
#        lattice = self.ek_lattice_class(
#            n_ghost_layers=1, agrid=self.params["agrid"])
#        ekspecies = self.make_default_ek_species(lattice)
#        # access out of bounds
#        out_of_bounds = max(ekspecies.shape) + 1
#        error_msg = 'Index error'
#        with self.assertRaisesRegex(Exception, error_msg):
#            ekspecies[out_of_bounds, 0, 0].density
#        with self.assertRaisesRegex(Exception, error_msg):
#            ekspecies[0, out_of_bounds, 0].density
#        with self.assertRaisesRegex(Exception, error_msg):
#            ekspecies[0, 0, out_of_bounds].density
#        # node index
#        node = ekspecies[1, 2, 3]
#        with self.assertRaisesRegex(RuntimeError, "Parameter 'index' is read-only"):
#            node.index = [2, 4, 6]
#        np.testing.assert_array_equal(np.copy(node.index), [1, 2, 3])
#        retval = node.call_method('override_index', index=[2, 4, 6])
#        self.assertEqual(retval, 0)
#        np.testing.assert_array_equal(np.copy(node.index), [2, 4, 6])
#        retval = node.call_method(
#            'override_index', index=[0, 0, out_of_bounds])
#        self.assertEqual(retval, 1)
#        np.testing.assert_array_equal(np.copy(node.index), [2, 4, 6])

    # TODO walberla: implement equality operator
#    def test_bool_operations_on_node(self):
#        lattice = self.ek_lattice_class(
#            n_ghost_layers=1, agrid=self.params["agrid"])
#        ekspecies = self.make_default_ek_species(lattice)
#        # test __eq()__ where a node is equal to itself and not equal to any
#        # other node
#        assert ekspecies[0, 0, 0] == ekspecies[0, 0, 0]
#        x, y, z = range(int(self.system.box_l[0])), range(
#            int(self.system.box_l[1])), range(int(self.system.box_l[2]))
#        nodes = [ekspecies[i, j, k] for i, j, k in itertools.product(x, y, z)]
#        nodes.remove(ekspecies[0, 0, 0])
#        assert all(ekspecies[0, 0, 0] != node for node in nodes)
#        # test __hash()__ intercept to identify nodes based on index rather
#        # than name. set() constructor runs hash()
#        subset1, subset2 = nodes[:-10], nodes[-10:]
#        assert len(set(subset1 + subset1)) == len(subset1)
#        assert len(set(subset1 + subset2)) == len(subset1) + len(subset2)


@utx.skipIfMissingFeatures("LB_WALBERLA")
class EKTestWalberla(EKTest, ut.TestCase):

    """Test for the Walberla implementation of the EK in double-precision."""

    ek_lattice_class = espressomd.lb.LatticeWalberla
    ek_species_class = espressomd.EKSpecies.EKSpecies
    ek_params = {"single_precision": False}
    atol = 1e-10
    rtol = 1e-7


@utx.skipIfMissingFeatures("LB_WALBERLA")
class EKTestWalberlaSinglePrecision(EKTest, ut.TestCase):

    """Test for the Walberla implementation of the EK in single-precision."""

    ek_lattice_class = espressomd.lb.LatticeWalberla
    ek_species_class = espressomd.EKSpecies.EKSpecies
    ek_params = {"single_precision": True}
    atol = 1e-7
    rtol = 5e-5


if __name__ == "__main__":
    ut.main()
