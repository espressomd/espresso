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

import numpy as np
import unittest as ut
import unittest_decorators as utx

import espressomd
import espressomd.lb
import espressomd.shapes
import espressomd.electrokinetics


class EKBoundariesBase:
    system = espressomd.System(box_l=[10.0, 5.0, 5.0])
    system.cell_system.skin = 0.1
    ek_species_params = {"kT": 1.5,
                         "density": 0.85,
                         "valency": 0.0,
                         "diffusion": 0.1,
                         "advection": False,
                         "friction_coupling": False,
                         "tau": 1.0}

    wall_shape1 = espressomd.shapes.Wall(normal=[1., 0., 0.], dist=2.5)
    wall_shape2 = espressomd.shapes.Wall(normal=[-1., 0., 0.], dist=-7.5)

    def setUp(self):
        self.lattice = self.ek_lattice_class(agrid=0.5, n_ghost_layers=1)

    def tearDown(self):
        self.system.ekcontainer.clear()

    def make_default_ek_species(self):
        return self.ek_species_class(
            lattice=self.lattice,
            single_precision=self.ek_params["single_precision"],
            **self.ek_species_params)

    def check_boundary_flags(self, ek_species, attr, value1, value2):
        def generator(value, shape):
            value_grid = np.tile(value, shape)
            if value_grid.shape[-1] == 1:
                value_grid = np.squeeze(value_grid, axis=-1)
            return value_grid

        accessor = np.vectorize(
            lambda obj: np.copy(getattr(obj, attr)),
            signature=f"()->({'n' if attr == 'flux' else ''})")

        slice1 = ek_species[:5, :, :]
        slice2 = ek_species[15:, :, :]
        slice3 = ek_species[5:15, :, :]
        np.testing.assert_equal(np.copy(slice1.is_boundary), True)
        np.testing.assert_equal(np.copy(slice2.is_boundary), True)
        np.testing.assert_equal(np.copy(slice3.is_boundary), False)
        field = f"{attr}_boundary"

        np.testing.assert_allclose(accessor(np.copy(getattr(slice1, field))),
                                   generator(value1, [5, 10, 10, 1]))
        np.testing.assert_allclose(accessor(np.copy(getattr(slice2, field))),
                                   generator(value2, [5, 10, 10, 1]))
        getattr(ek_species, f"clear_{attr}_boundaries")()
        np.testing.assert_equal(
            np.copy(ek_species[:, :, :].is_boundary), False)

    def test_flux_boundary_flags(self):
        flux1 = 1e-3 * np.array([1., 2., 3.])
        flux2 = 1e-3 * np.array([4., 5., 6.])

        # check with two shapes
        ek_species = self.make_default_ek_species()
        value_shape = tuple(ek_species.shape) + (3,)
        ek_species.add_boundary_from_shape(
            shape=self.wall_shape1, value=flux1,
            boundary_type=espressomd.electrokinetics.FluxBoundary)
        ek_species.add_boundary_from_shape(
            shape=self.wall_shape2, value=flux2 * np.ones(value_shape),
            boundary_type=espressomd.electrokinetics.FluxBoundary)
        self.check_boundary_flags(ek_species, "flux", flux1, flux2)

        # check with union of two shapes
        ek_species = self.make_default_ek_species()
        union = espressomd.shapes.Union()
        union.add([self.wall_shape1, self.wall_shape2])
        ek_species.add_boundary_from_shape(
            shape=union, value=flux1,
            boundary_type=espressomd.electrokinetics.FluxBoundary)
        self.check_boundary_flags(ek_species, "flux", flux1, flux1)

    def test_density_boundary_flags(self):
        density1 = 1.
        density2 = 2.

        # check with two shapes
        ek_species = self.make_default_ek_species()
        value_shape = tuple(ek_species.shape) + (1,)
        ek_species.add_boundary_from_shape(
            shape=self.wall_shape1, value=density1,
            boundary_type=espressomd.electrokinetics.DensityBoundary)
        ek_species.add_boundary_from_shape(
            shape=self.wall_shape2, value=density2 * np.ones(value_shape),
            boundary_type=espressomd.electrokinetics.DensityBoundary)
        self.check_boundary_flags(ek_species, "density", density1, density2)

        # check with union of two shapes
        ek_species = self.make_default_ek_species()
        union = espressomd.shapes.Union()
        union.add([self.wall_shape1, self.wall_shape2])
        ek_species.add_boundary_from_shape(
            shape=union, value=density1,
            boundary_type=espressomd.electrokinetics.DensityBoundary)
        self.check_boundary_flags(ek_species, "density", density1, density1)

    def test_exceptions(self):
        ek_species = self.make_default_ek_species()
        with self.assertRaisesRegex(TypeError, "Parameter 'boundary_type' must be a subclass of FluxBoundary or DensityBoundary"):
            ek_species.add_boundary_from_shape(
                shape=self.wall_shape1, value=[0., 0., 0.],
                boundary_type=espressomd.lb.VelocityBounceBack)
        with self.assertRaisesRegex(ValueError, "expected an espressomd.shapes.Shape"):
            ek_species.add_boundary_from_shape(
                shape=ek_species, value=[0., 0., 0.],
                boundary_type=espressomd.electrokinetics.FluxBoundary)
        with self.assertRaisesRegex(ValueError, r"Cannot process density value grid of shape \(3,\)"):
            ek_species.add_boundary_from_shape(
                shape=self.wall_shape1, value=[0., 0., 0.],
                boundary_type=espressomd.electrokinetics.DensityBoundary)
        with self.assertRaisesRegex(ValueError, r"Cannot process flux value grid of shape \(1,\)"):
            ek_species.add_boundary_from_shape(
                shape=self.wall_shape1, value=0.,
                boundary_type=espressomd.electrokinetics.FluxBoundary)


@utx.skipIfMissingFeatures(["WALBERLA"])
class EKBoundariesWalberla(EKBoundariesBase, ut.TestCase):

    """Test for the Walberla implementation of the LB in double-precision."""

    ek_lattice_class = espressomd.electrokinetics.LatticeWalberla
    ek_species_class = espressomd.electrokinetics.EKSpecies
    ek_params = {"single_precision": False}


@utx.skipIfMissingFeatures(["WALBERLA"])
class EKBoundariesWalberlaSinglePrecision(EKBoundariesBase, ut.TestCase):

    """Test for the Walberla implementation of the LB in single-precision."""

    ek_lattice_class = espressomd.electrokinetics.LatticeWalberla
    ek_species_class = espressomd.electrokinetics.EKSpecies
    ek_params = {"single_precision": True}


if __name__ == "__main__":
    ut.main()
