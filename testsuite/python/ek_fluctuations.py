#
# Copyright (C) 2024-2026 The ESPResSo project
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
import itertools


class EKTest:
    TAU = 1.123
    DENSITY = 27.0
    DIFFUSION_COEFFICIENT = 0.678
    AGRID = 2.865
    BOX_L = 8 * AGRID

    system = espressomd.System(box_l=[BOX_L, BOX_L, BOX_L])
    system.time_step = TAU
    system.cell_system.skin = 0.4

    def tearDown(self):
        self.system.ekcontainer.clear()

    def test_fluctuation(self):
        decimal_precision = 1 if self.lattice_params["single_precision"] else 9

        lattice = espressomd.electrokinetics.Lattice(
            n_ghost_layers=1, agrid=self.AGRID)

        target_particlenum = self.DENSITY * self.system.volume()

        species = espressomd.electrokinetics.EKSpecies(
            lattice=lattice, density=self.DENSITY, valency=0.0, advection=False,
            diffusion=self.DIFFUSION_COEFFICIENT, friction_coupling=False,
            tau=self.TAU, thermalized=True, seed=42, **self.lattice_params)

        eksolver = espressomd.electrokinetics.EKNone(
            lattice=lattice, tau=self.TAU)
        self.system.ekcontainer = espressomd.electrokinetics.EKContainer(
            tau=self.TAU, solver=eksolver)
        self.system.ekcontainer.add(species)

        self.system.integrator.run(100)

        n_min = 10.0
        n_max = 44.0
        bin_size = 0.25
        x_range = np.linspace(
            n_min + 0.5 * bin_size, n_max - 0.5 * bin_size,
            int((n_max - n_min) / bin_size))
        sample_steps = 150
        integration_steps = 75

        bins = int((n_max - n_min) / bin_size)
        hist, _ = np.histogram(
            [], bins=bins, range=(n_min, n_max), density=False)

        cell_volume = self.AGRID**3

        for _ in range(sample_steps):
            self.system.integrator.run(integration_steps)

            local_density = species[:, :, :].density

            np.testing.assert_almost_equal(
                np.sum(local_density * cell_volume), target_particlenum, decimal=decimal_precision)

            hist += np.histogram(local_density, bins=bins,
                                 range=(n_min, n_max), density=False)[0]

        hist = hist / np.sum(hist) / bin_size

        # Positive half of the D3Q27 neighbor set used by the staggered EK flux.
        directions = np.array(
            [offset for offset in itertools.product(
                (-1, 0, 1), repeat=3) if offset > (0, 0, 0)],
            dtype=int
        )
        A0 = sum([np.linalg.norm(offset) for offset in directions]) / 3.0
        weights = 1 / \
            (A0 * np.array([np.linalg.norm(offset) for offset in directions]))

        # get all fourier-vectors of the lattice, excluding the zero mode
        k_values = np.indices(lattice.shape).reshape(3, -1).T
        k_vectors = 2.0 * np.pi * k_values / lattice.shape
        k_vectors = k_vectors[np.any(k_values != 0, axis=1)]

        diffusion_lattice = self.DIFFUSION_COEFFICIENT * self.TAU / self.AGRID**2
        # eigen-value contribution per timestep of the discrete diffusion operator
        lambdas = -2.0 * diffusion_lattice * np.sum(
            weights[None, :] * (1.0 - np.cos(k_vectors @ directions.T)),
            axis=1
        )
        # explicit Euler causes amplification of the poisson-distributed noise
        # take average over amplification of all modes
        amplification = np.mean(1.0 / (1.0 + 0.5 * lambdas))
        effective_scale = cell_volume / amplification

        mean_particles = self.DENSITY * effective_scale
        particles_per_cell = x_range * effective_scale
        analytic_distribution = effective_scale * 1 / np.sqrt(2.0 * np.pi * particles_per_cell) * np.power(
            mean_particles / particles_per_cell, particles_per_cell) * np.exp(particles_per_cell - mean_particles)

        np.testing.assert_allclose(
            hist, analytic_distribution, rtol=0, atol=0.008)


@utx.skipIfMissingFeatures(["WALBERLA"])
class EKFluctuationsDoublePrecisionCPU(EKTest, ut.TestCase):
    lattice_params = {"single_precision": False, "gpu": False}


@utx.skipIfMissingFeatures(["WALBERLA"])
class EKFluctuationsSinglePrecisionCPU(EKTest, ut.TestCase):
    lattice_params = {"single_precision": True, "gpu": False}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class EKFluctuationsDoublePrecisionGPU(EKTest, ut.TestCase):
    lattice_params = {"single_precision": False, "gpu": True}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class EKFluctuationsSinglePrecisionGPU(EKTest, ut.TestCase):
    lattice_params = {"single_precision": True, "gpu": True}


if __name__ == "__main__":
    ut.main()
