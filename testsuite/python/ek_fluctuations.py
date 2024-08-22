#
# Copyright (C) 2024 The ESPResSo project
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
import math


@utx.skipIfMissingFeatures(["WALBERLA"])
class EKFluctuations(ut.TestCase):
    BOX_L = 8
    TAU = 0.1
    DENSITY = 27.0
    DIFFUSION_COEFFICIENT = 0.01
    AGRID = 1.0

    system = espressomd.System(box_l=[BOX_L, BOX_L, BOX_L])
    system.time_step = TAU
    system.cell_system.skin = 0.4

    def tearDown(self) -> None:
        self.system.ekcontainer.clear()

    def test_diffusion_single(self):
        self.detail_test_fluctuation(single_precision=True)

    def test_diffusion_double(self):
        self.detail_test_fluctuation(single_precision=False)

    def detail_test_fluctuation(self, single_precision: bool):
        decimal_precision: int = 2 if single_precision else 10

        lattice = espressomd.electrokinetics.LatticeWalberla(
            n_ghost_layers=1, agrid=self.AGRID)

        target_density = self.DENSITY * self.system.volume()

        species = espressomd.electrokinetics.EKSpecies(
            lattice=lattice, density=self.DENSITY, valency=0.0, advection=False,
            diffusion=self.DIFFUSION_COEFFICIENT, friction_coupling=False,
            single_precision=single_precision, tau=self.TAU, thermalized=True, seed=42)

        eksolver = espressomd.electrokinetics.EKNone(lattice=lattice)

        self.system.ekcontainer = espressomd.electrokinetics.EKContainer(
            tau=self.TAU, solver=eksolver)
        self.system.ekcontainer.add(species)

        self.system.integrator.run(100)

        # Set integration and binning parameters
        n_min = 10.0
        n_max = 44.0
        bin_size = 0.25
        x_range = np.linspace(n_min, n_max, int((n_max - n_min) / bin_size))
        sample_steps = 150
        integration_steps = 100

        bins = int((n_max - n_min) / bin_size)
        hist, _ = np.histogram(
            [], bins=bins, range=(n_min, n_max), density=False)

        # Integrate
        for _ in range(sample_steps):
            self.system.integrator.run(integration_steps)

            density = species[:, :, :].density

            np.testing.assert_almost_equal(
                np.sum(density), target_density, decimal=decimal_precision)

            hist += np.histogram(density, bins=bins,
                                 range=(n_min, n_max), density=False)[0]

        hist = hist / np.sum(hist) / bin_size

        analytic_distribution = 1.0 / np.sqrt(2.0 * math.pi * x_range) * np.power(
            self.DENSITY / x_range, x_range) * np.exp(x_range - self.DENSITY)

        max_diff = np.max(np.abs(analytic_distribution - hist))
        self.assertLess(max_diff, 5.0e-03,
                        f"Density distribution accuracy not achieved, allowed "
                        f"deviation: 5.0e-03, measured: {max_diff}")


if __name__ == "__main__":
    ut.main()
