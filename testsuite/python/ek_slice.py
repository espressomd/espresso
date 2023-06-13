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

import itertools
import numpy as np
import unittest as ut
import unittest_decorators as utx

import espressomd
import espressomd.electrokinetics


@utx.skipIfMissingFeatures("WALBERLA")
class Test(ut.TestCase):

    """This simple test first writes random numbers and then reads them
    to same slices of LB nodes and compares if the results are the same,
    shape-wise and value-wise.
    """

    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    system.time_step = .01
    system.cell_system.skin = 0.1
    ek_species_params = {"kT": 1.5,
                         "density": 0.85,
                         "valency": 0.0,
                         "diffusion": 0.1,
                         "advection": False,
                         "friction_coupling": False,
                         "tau": 1.0}
    np.random.seed(seed=42)

    @classmethod
    def setUpClass(cls):
        cls.lattice = espressomd.electrokinetics.LatticeWalberla(
            agrid=1., n_ghost_layers=1)
        cls.ek_species = espressomd.electrokinetics.EKSpecies(
            lattice=cls.lattice,
            single_precision=False,
            **cls.ek_species_params)

    def test_slicing(self):
        ek_species = self.ek_species

        # array locked
        array = ek_species[1:-1:1, 5, 3:6].density
        with self.assertRaisesRegex(ValueError, "ESPResSo array properties return non-writable arrays"):
            array[0, 0, 0] = 5.

        input_dens = np.random.rand(8, 3) + 1.
        ek_species[1:-1, 5, 3:6].density = input_dens
        output_dens = ek_species[1:-1, 5, 3:6].density
        np.testing.assert_array_almost_equal(np.copy(output_dens), input_dens)

        # density broadcast (with type conversion from int to double)
        ek_species[:, :, 0].density = 2
        np.testing.assert_array_almost_equal(
            np.copy(ek_species[:, :, 0].density), 2.)

        # flux boundary on slice
        output_boundary_shape = ek_species[1:, 1:, 1:].flux_boundary.shape
        should_boundary_shape = (9, 9, 9)
        np.testing.assert_array_equal(
            output_boundary_shape, should_boundary_shape)

        with self.assertRaisesRegex(TypeError, "Parameter 'values' must be an array_like of FluxBoundary or None"):
            ek_species[1:, 1:, 1:].flux_boundary = np.zeros(
                should_boundary_shape)
        with self.assertRaisesRegex(TypeError, "Parameter 'values' must be an array_like of FluxBoundary or None"):
            ek_species[1:, 1:, 1:].flux_boundary = np.array(
                [None, [1, 2, 3]], dtype=object)

        flux_ref = espressomd.electrokinetics.FluxBoundary([1e-6, 2e-6, 3e-6])
        ek_species[1:2, 1:, 0].flux_boundary = flux_ref
        ek_species[1:2, 2:, 0].flux_boundary = None
        for flux in ek_species[1:2, 1, 0].flux_boundary.flatten():
            np.testing.assert_array_almost_equal(
                flux.flux, flux_ref.flux)
        for flux in ek_species[1:2, 2:, 0:2].flux_boundary.flatten():
            self.assertIsNone(flux)

        # density boundary on slice
        output_boundary_shape = ek_species[0:-1, 1:, 1:].density_boundary.shape
        should_boundary_shape = (9, 9, 9)
        np.testing.assert_array_equal(
            output_boundary_shape, should_boundary_shape)

        with self.assertRaisesRegex(TypeError, "Parameter 'values' must be an array_like of DensityBoundary or None"):
            ek_species[1:, 1:, 1:].density_boundary = np.zeros(
                should_boundary_shape)
        with self.assertRaisesRegex(TypeError, "Parameter 'values' must be an array_like of DensityBoundary or None"):
            ek_species[1:, 1:, 1:].density_boundary = np.array(
                [None, 1.], dtype=object)

        dens_ref = espressomd.electrokinetics.DensityBoundary(1e-6)
        ek_species[2:3, 1:, 0].density_boundary = dens_ref
        ek_species[2:3, 2:, 0].density_boundary = None
        for dens in ek_species[2:3, 1, 0].density_boundary.flatten():
            np.testing.assert_array_almost_equal(
                dens.density, dens_ref.density)
        for dens in ek_species[2:3, 2:, 0:2].density_boundary.flatten():
            self.assertIsNone(dens)

        # is_boundary on slice
        output_boundary_shape = ek_species[1:, 1:, 1:].is_boundary.shape
        should_boundary_shape = (9, 9, 9)
        np.testing.assert_array_equal(
            output_boundary_shape, should_boundary_shape)
        np.testing.assert_array_equal(
            np.copy(ek_species[1:3, 1:2, 0].is_boundary), True)
        np.testing.assert_array_equal(
            np.copy(ek_species[3:, 2:, 0:2].is_boundary), False)

        with self.assertRaisesRegex(RuntimeError, "Property 'is_boundary' is read-only"):
            ek_species[1:, 1:, 1:].is_boundary = np.zeros(
                should_boundary_shape)

        # access out of bounds
        i = ek_species.shape[2] + 10
        ek_slice = ek_species[1, 2, i:i + 10]
        self.assertEqual(ek_slice.density.shape, (0,))
        self.assertIsInstance(ek_slice.density.dtype, object)
        with self.assertRaisesRegex(AttributeError, "Cannot set properties of an empty 'EKSpeciesSlice' object"):
            ek_slice.density = [1., 2., 3.]

        # other exceptions
        with self.assertRaisesRegex(RuntimeError, "Unknown EK property 'unknown'"):
            ek_species[:, :, :].call_method("get_value_shape", name="unknown")

    def test_iterator(self):
        ekslice_handle = self.ek_species[:, :, :]
        # arrange node indices using class methods
        ek_indices = [np.arange(self.ek_species.shape[i]) for i in range(3)]
        arranged_indices = list(itertools.product(*ek_indices))
        # arrange node indices using __iter__() enforced conversion
        iterator_indices = [node.index for node in ekslice_handle]
        # check the results correspond pairwise. order is implicitly preserved.
        np.testing.assert_array_equal(arranged_indices, iterator_indices)
        # use __eq()__ method form EKSpeciesNode()
        assert all([x == y for x, y in zip(
            arranged_indices, iterator_indices)])


if __name__ == "__main__":
    ut.main()
