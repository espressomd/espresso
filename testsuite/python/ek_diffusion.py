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

import unittest as ut
import unittest_decorators as utx
import espressomd
import espressomd.electrokinetics
import numpy as np
import scipy.optimize


@utx.skipIfMissingFeatures(["WALBERLA"])
class EKDiffusion(ut.TestCase):
    BOX_L = 15.5
    AGRID = 0.5
    DENSITY = 1
    DIFFUSION_COEFFICIENT = 0.05
    TAU = 0.9
    TIMESTEPS = int(65 / TAU)

    system = espressomd.System(box_l=[BOX_L, BOX_L, BOX_L])
    system.time_step = TAU
    system.cell_system.skin = 0.4

    def tearDown(self) -> None:
        self.system.ekcontainer.clear()

    def analytical_density(self, pos: np.ndarray, time: int, D: float):
        return (4 * np.pi * D * time)**(-3 / 2) * \
            np.exp(-np.sum(np.square(pos), axis=-1) / (4 * D * time))

    def test_diffusion_single(self):
        self.detail_test_diffusion(single_precision=True)

    def test_diffusion_double(self):
        self.detail_test_diffusion(single_precision=False)

    def detail_test_diffusion(self, single_precision: bool):
        """
        Testing EK for simple diffusion of a point droplet
        """

        decimal_precision: int = 7 if single_precision else 10

        lattice = espressomd.electrokinetics.LatticeWalberla(
            n_ghost_layers=1, agrid=self.AGRID)

        ekspecies = espressomd.electrokinetics.EKSpecies(
            lattice=lattice, density=0.0, valency=0.0, advection=False,
            diffusion=self.DIFFUSION_COEFFICIENT, friction_coupling=False,
            single_precision=single_precision, tau=self.TAU)

        eksolver = espressomd.electrokinetics.EKNone(lattice=lattice)

        self.system.ekcontainer.tau = self.TAU
        self.system.ekcontainer.solver = eksolver
        self.system.ekcontainer.add(ekspecies)

        center = np.asarray(lattice.shape // 2, dtype=int)

        ekspecies[center].density = self.DENSITY

        # check that the density in the domain is what is expected
        np.testing.assert_almost_equal(
            np.sum(ekspecies[:, :, :].density), self.DENSITY, decimal_precision)

        # calculate physical positions
        positions = np.empty((*lattice.shape, 3))
        positions[..., 2], positions[..., 1], positions[..., 0] = np.meshgrid(
            *map(lambda x: np.arange(0, x) - x / 2, lattice.shape))
        positions += 0.5
        positions *= self.AGRID

        self.system.integrator.run(self.TIMESTEPS)

        simulated_density = np.copy(ekspecies[:, :, :].density)

        # check that the density is conserved
        np.testing.assert_almost_equal(
            np.sum(simulated_density), self.DENSITY, decimal_precision)
        assert np.all(simulated_density >= 0.), "EK has negative densities"

        # check that the maximum is in the right place
        peak = np.unravel_index(
            np.argmax(simulated_density, axis=None),
            lattice.shape)
        np.testing.assert_equal(peak, center)

        calc_density = self.analytical_density(
            positions, self.TIMESTEPS * self.TAU, self.DIFFUSION_COEFFICIENT) * self.AGRID ** 3

        target = [self.TIMESTEPS * self.TAU, self.DIFFUSION_COEFFICIENT]

        popt, _ = scipy.optimize.curve_fit(self.analytical_density,
                                           positions.reshape(-1, 3),
                                           simulated_density.reshape(
                                               -1) / self.AGRID ** 3,
                                           p0=target,
                                           bounds=([0, 0], [np.inf, np.inf]))

        np.testing.assert_allclose(
            popt[0], self.TIMESTEPS * self.TAU, rtol=0.1)
        np.testing.assert_allclose(
            popt[1], self.DIFFUSION_COEFFICIENT, rtol=0.1)
        np.testing.assert_allclose(
            calc_density, simulated_density, atol=1e-5, rtol=0.)


if __name__ == "__main__":
    ut.main()
