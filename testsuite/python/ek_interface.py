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

    def tearDown(self):
        self.system.ekcontainer.clear()
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.time_step = self.params["tau"]

    def make_default_ek_species(self, lattice):
        return self.ek_species_class(
            lattice=lattice,
            single_precision=self.ek_params["single_precision"],
            **self.ek_species_params)

    def test_properties(self):
        lattice = self.ek_lattice_class(
            n_ghost_layers=1, agrid=self.params["agrid"])

        # inactive actor
        ekspecies = self.make_default_ek_species(lattice)
        self.check_properties(ekspecies)

        eksolver = espressomd.EKSpecies.EKNone(lattice=lattice)
        self.system.ekcontainer.tau = self.system.time_step
        self.system.ekcontainer.solver = eksolver
        self.assertAlmostEqual(
            self.system.ekcontainer.tau,
            self.system.time_step,
            delta=self.atol)

        # activated actor
        ekspecies = self.make_default_ek_species(lattice)
        self.system.ekcontainer.add(ekspecies)
        self.check_properties(ekspecies)

        # deactivated actor
        ekspecies = self.make_default_ek_species(lattice)
        self.system.ekcontainer.add(ekspecies)
        self.system.ekcontainer.remove(ekspecies)
        self.check_properties(ekspecies)

    def check_properties(self, species):
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

    def test_raise_if_read_only(self):
        lattice = self.ek_lattice_class(
            n_ghost_layers=1, agrid=self.params["agrid"])
        ekspecies = self.make_default_ek_species(lattice)
        for key in {"lattice", "shape", "is_single_precision"}:
            with self.assertRaisesRegex(RuntimeError, f"(Parameter|Property) '{key}' is read-only"):
                setattr(ekspecies, key, 0)

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
