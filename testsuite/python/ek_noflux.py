#
# Copyright (C) 2022-2026 The ESPResSo project
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
import espressomd.shapes
import espressomd.electrokinetics


class EKTest:
    AGRID = 2.865
    BOX_L = 16. * AGRID
    DENSITY = 1
    DIFFUSION_COEFFICIENT = 0.1
    TIME = 50
    RADIUS = 5.

    system = espressomd.System(box_l=[BOX_L, BOX_L, BOX_L])
    TAU = 1.85
    system.time_step = TAU
    system.cell_system.skin = 0.4

    def tearDown(self):
        self.system.ekcontainer = None

    def test_noflux(self):
        """
        Testing the EK noflux boundaries to not leak density outside of a sphere.
        """

        decimal_precision = 7 if self.lattice_params["single_precision"] else 10

        lattice = espressomd.electrokinetics.Lattice(
            n_ghost_layers=2, agrid=self.AGRID)

        ekspecies = espressomd.electrokinetics.EKSpecies(
            lattice=lattice, density=0.0, diffusion=self.DIFFUSION_COEFFICIENT,
            valency=0.0, advection=False, friction_coupling=False,
            tau=self.TAU, **self.lattice_params)

        eksolver = espressomd.electrokinetics.EKNone(
            lattice=lattice, tau=self.TAU)

        self.system.ekcontainer = espressomd.electrokinetics.EKContainer(
            tau=self.TAU, solver=eksolver)
        self.system.ekcontainer.add(ekspecies)

        center = np.asarray(lattice.shape / 2, dtype=int)

        ekspecies[center[0], center[1], center[2]].density = self.DENSITY

        sphere = espressomd.shapes.Sphere(
            center=self.system.box_l / 2,
            radius=self.RADIUS,
            direction=-1)
        ekspecies.add_boundary_from_shape(
            sphere, [0, 0, 0], espressomd.electrokinetics.FluxBoundary)

        positions = np.empty((*lattice.shape, 3))
        positions[..., 2], positions[..., 1], positions[..., 0] = np.meshgrid(
            *map(lambda x: np.arange(0, x) - x / 2, lattice.shape))
        positions += 0.5
        positions *= self.AGRID

        self.system.integrator.run(self.TIME)

        simulated_density = np.copy(ekspecies[:, :, :].density)

        # check that the density is conserved globally
        np.testing.assert_almost_equal(
            np.sum(simulated_density), self.DENSITY, decimal_precision)

        domain_density = simulated_density[np.logical_not(
            ekspecies[:, :, :].is_boundary)]
        # check that the density is kept constant inside the sphere
        np.testing.assert_almost_equal(
            np.sum(domain_density), self.DENSITY, decimal_precision)
        np.testing.assert_array_less(
            0., domain_density, "EK density array contains negative densities!")


@utx.skipIfMissingFeatures(["WALBERLA"])
class EKNoFluxDoublePrecisionCPU(EKTest, ut.TestCase):
    lattice_params = {"single_precision": False, "gpu": False}


@utx.skipIfMissingFeatures(["WALBERLA"])
class EKNoFluxSinglePrecisionCPU(EKTest, ut.TestCase):
    lattice_params = {"single_precision": True, "gpu": False}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class EKNoFluxDoublePrecisionGPU(EKTest, ut.TestCase):
    lattice_params = {"single_precision": False, "gpu": True}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class EKNoFluxSinglePrecisionGPU(EKTest, ut.TestCase):
    lattice_params = {"single_precision": True, "gpu": True}


if __name__ == "__main__":
    ut.main()
