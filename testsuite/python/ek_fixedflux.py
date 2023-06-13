#
# Copyright (C) 2022-2023 The ESPResSo project
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


@utx.skipIfMissingFeatures(["WALBERLA"])
class EKFixedFlux(ut.TestCase):
    BOX_L = 2.5
    AGRID = 0.5
    DENSITY = 1
    DIFFUSION_COEFFICIENT = 0.1
    TIMESTEPS = 40
    TAU = 0.25

    INFLOW_FLUX = 0.1

    system = espressomd.System(box_l=[BOX_L, BOX_L, BOX_L])
    system.time_step = 1.0
    system.cell_system.skin = 0.4

    def tearDown(self) -> None:
        self.system.ekcontainer.clear()

    def test_inflow_single(self):
        self.detail_test_inflow(single_precision=True)

    def test_inflow_double(self):
        self.detail_test_inflow(single_precision=False)

    def detail_test_inflow(self, single_precision: bool):
        """
        Testing the EK fixed flux boundaries to test the fixed inflow into a non-periodic box.
        """

        decimal_precision: int = 5 if single_precision else 7

        lattice = espressomd.electrokinetics.LatticeWalberla(
            n_ghost_layers=1, agrid=self.AGRID)

        ekspecies = espressomd.electrokinetics.EKSpecies(
            lattice=lattice, density=0.0, diffusion=self.DIFFUSION_COEFFICIENT,
            valency=0.0, advection=False, friction_coupling=False,
            single_precision=single_precision, tau=self.TAU)

        eksolver = espressomd.electrokinetics.EKNone(lattice=lattice)

        self.system.ekcontainer.tau = self.TAU
        self.system.ekcontainer.solver = eksolver
        self.system.ekcontainer.add(ekspecies)

        ekspecies[1:-1, 1:-1, 1:-1].density = self.DENSITY

        ekspecies[:, :, 0].flux_boundary = \
            espressomd.electrokinetics.FluxBoundary([0, 0, 0])
        ekspecies[:, :, -1].flux_boundary = \
            espressomd.electrokinetics.FluxBoundary([0, 0, 0])
        ekspecies[:, 0, :].flux_boundary = \
            espressomd.electrokinetics.FluxBoundary([0, 0, 0])
        ekspecies[:, -1, :].flux_boundary = \
            espressomd.electrokinetics.FluxBoundary([0, 0, 0])
        ekspecies[0, :, :].flux_boundary = \
            espressomd.electrokinetics.FluxBoundary([0, 0, 0])
        ekspecies[-1, :, :].flux_boundary = \
            espressomd.electrokinetics.FluxBoundary([0, 0, 0])

        # set fixed flux in +z-direction
        ekspecies[:, :, -1].flux_boundary = espressomd.electrokinetics.FluxBoundary(
            [0, 0, -self.INFLOW_FLUX])
        additional_center_flux = 3 * self.INFLOW_FLUX
        midpoint = int(lattice.shape[0] // 2)
        ekspecies[midpoint, midpoint, -1].flux_boundary = \
            espressomd.electrokinetics.FluxBoundary(
                [0, 0, -self.INFLOW_FLUX - additional_center_flux])

        # check density before integration
        expected_initial_density = self.DENSITY * np.prod(lattice.shape - 2)

        np.testing.assert_almost_equal(
            actual=np.sum(ekspecies[1:-1, 1:-1, 1:-1].density),
            desired=expected_initial_density, decimal=decimal_precision)

        self.system.integrator.run(self.TIMESTEPS)

        # check that density has pushed into domain
        inflow_area = np.prod(lattice.shape[:2] - 2)
        expected_end_density = expected_initial_density + \
            (self.INFLOW_FLUX * inflow_area + additional_center_flux) * \
            self.TIMESTEPS * self.TAU / self.AGRID

        np.testing.assert_almost_equal(
            actual=np.sum(ekspecies[1:-1, 1:-1, 1:-1].density),
            desired=expected_end_density, decimal=decimal_precision)


if __name__ == "__main__":
    ut.main()
