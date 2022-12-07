#
# Copyright (C) 2022 The ESPResSo project
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
import espressomd.EKSpecies
import numpy as np
import scipy.optimize


@utx.skipIfMissingFeatures(["WALBERLA"])
class EKDiffusion(ut.TestCase):
    BOX_L = 31.
    AGRID = 1.0
    DENSITY = 1
    DIFFUSION_COEFFICIENT = 0.1
    TIME = 150

    system = espressomd.System(box_l=[BOX_L, BOX_L, BOX_L])
    system.time_step = 1.0
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

        lattice = espressomd.lb.LatticeWalberla(
            n_ghost_layers=1, agrid=self.AGRID)

        ekspecies = espressomd.EKSpecies.EKSpecies(lattice=lattice,
                                                   density=0.0, kT=0.0, diffusion=self.DIFFUSION_COEFFICIENT, valency=0.0,
                                                   advection=False, friction_coupling=False, ext_efield=[0, 0, 0], single_precision=single_precision)

        eksolver = espressomd.EKSpecies.EKNone(lattice=lattice)

        self.system.ekcontainer.add(ekspecies)
        self.system.ekcontainer.tau = 1.0
        self.system.ekcontainer.solver = eksolver

        center = np.asarray(self.system.box_l / 2, dtype=np.int)

        ekspecies[center].density = self.DENSITY

        # check that the density in the domain is what is expected
        np.testing.assert_almost_equal(
            np.sum(ekspecies[:, :, :].density), self.DENSITY, decimal_precision)

        # TODO: replace that when the blockforest is able to return the
        # dimensions
        positions = np.empty((*self.system.box_l.astype(np.int), 3))
        positions[..., 2], positions[..., 1], positions[..., 0] = np.meshgrid(
            *map(lambda x: np.arange(0, x) - x / 2, self.system.box_l))
        positions += 0.5

        self.system.integrator.run(self.TIME)

        simulated_density = np.copy(ekspecies[:, :, :].density)

        # check that the density is conserved
        np.testing.assert_almost_equal(
            np.sum(simulated_density), self.DENSITY, decimal_precision)
        if np.any(simulated_density < 0.):
            self.fail("EK density array contains negative densities!")

        # check that the maximum is in the right place
        peak = np.unravel_index(
            np.argmax(simulated_density, axis=None),
            self.system.box_l.astype(np.int))
        np.testing.assert_equal(peak, self.system.box_l / 2 - 0.5)

        calc_density = self.analytical_density(
            positions, self.TIME, self.DIFFUSION_COEFFICIENT)
        target = [self.TIME, self.DIFFUSION_COEFFICIENT]

        popt, _ = scipy.optimize.curve_fit(self.analytical_density,
                                           positions.reshape(-1, 3),
                                           simulated_density.reshape(-1),
                                           p0=target,
                                           bounds=([0, 0], [np.inf, np.inf]))

        np.testing.assert_allclose(popt[0], self.TIME, rtol=0.1)
        np.testing.assert_allclose(
            popt[1], self.DIFFUSION_COEFFICIENT, rtol=0.1)
        np.testing.assert_allclose(
            calc_density,
            simulated_density,
            atol=1e-5,
            rtol=0)


if __name__ == "__main__":
    ut.main()
