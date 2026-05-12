#
# Copyright (C) 2026 The ESPResSo project
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


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class FFTPotential:
    BOX_L = [22.92, 28.65, 34.38]
    AGRID = 2.865
    TAU = 5.943635

    system = espressomd.System(box_l=BOX_L)
    system.time_step = TAU
    system.cell_system.skin = 0.4

    def tearDown(self):
        self.system.ekcontainer = None

    def test_potential(self):
        eps0 = 0.09135
        epsR = 8.0
        kT = 3.15379
        valency = 1.9
        rho0 = 1.3

        external_electric_field = np.asarray([0.0, 0.0, 0.0])

        lattice = espressomd.electrokinetics.Lattice(
            n_ghost_layers=1, agrid=self.AGRID)

        ekspecies_pos = espressomd.electrokinetics.EKSpecies(
            lattice=lattice, density=0.0, kT=kT, valency=valency,
            diffusion=0.0, friction_coupling=False,
            advection=False, ext_efield=external_electric_field,
            tau=self.TAU, **self.lattice_params)

        ekspecies_neg = espressomd.electrokinetics.EKSpecies(
            lattice=lattice, density=0.0, kT=kT, valency=-valency,
            diffusion=0.0, friction_coupling=False,
            advection=False, ext_efield=external_electric_field,
            tau=self.TAU, **self.lattice_params)

        eksolver = espressomd.electrokinetics.EKFFT(
            lattice=lattice, permittivity=eps0 * epsR, tau=self.TAU, **self.lattice_params)

        self.system.ekcontainer = espressomd.electrokinetics.EKContainer(
            tau=self.TAU, solver=eksolver)
        self.system.ekcontainer.add(ekspecies_pos)
        self.system.ekcontainer.add(ekspecies_neg)

        Nx, Ny, Nz = lattice.shape
        Lx, Ly, Lz = self.system.box_l
        nx, ny, nz = (1, 2, 4)

        hx = Lx / Nx
        hy = Ly / Ny
        hz = Lz / Nz

        i = np.arange(Nx, dtype=float)
        j = np.arange(Ny, dtype=float)
        k = np.arange(Nz, dtype=float)

        I, J, K = np.meshgrid(i, j, k, indexing="ij")
        phase = 2.0 * np.pi * (nx * I / Nx + ny * J / Ny + nz * K / Nz)

        phi_exact = rho0 * np.cos(phase)

        lap = (
            (np.roll(phi_exact, -1, axis=0) - 2.0 * phi_exact + np.roll(phi_exact, 1, axis=0)) / hx**2 +
            (np.roll(phi_exact, -1, axis=1) - 2.0 * phi_exact + np.roll(phi_exact, 1, axis=1)) / hy**2 +
            (np.roll(phi_exact, -1, axis=2) - 2.0 * phi_exact + np.roll(phi_exact, 1, axis=2)) / hz**2  # nopep8
        )

        rho = - eps0 * epsR / valency * lap
        rho -= rho.mean()
        phi_exact -= phi_exact.mean()

        dens_pos = np.zeros(lattice.shape, dtype=float)
        dens_neg = np.zeros(lattice.shape, dtype=float)
        for i, j, k in itertools.product(range(Nx), range(Ny), range(Nz)):
            dens_pos[i, j, k] = max(+rho[i, j, k], 0.0)
            dens_neg[i, j, k] = max(-rho[i, j, k], 0.0)
        ekspecies_pos[:, :, :].density = dens_pos
        ekspecies_neg[:, :, :].density = dens_neg
        self.system.integrator.run(1)

        calc_potential = np.copy(eksolver[:, :, :].potential)
        calc_potential -= calc_potential.mean()

        np.testing.assert_allclose(
            calc_potential, phi_exact, rtol=0.0, atol=self.atol)


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class EKTestWalberla(FFTPotential, ut.TestCase):
    lattice_params = {"single_precision": False, "gpu": False}
    atol = 5e-12


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class EKTestWalberlaSinglePrecision(FFTPotential, ut.TestCase):
    lattice_params = {"single_precision": True, "gpu": False}
    atol = 5e-7


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT", "CUDA"])
class EKTestWalberlaGPU(FFTPotential, ut.TestCase):
    lattice_params = {"single_precision": False, "gpu": True}
    atol = 5e-12


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT", "CUDA"])
class EKTestWalberlaSinglePrecisionGPU(FFTPotential, ut.TestCase):
    lattice_params = {"single_precision": True, "gpu": True}
    atol = 4e-7


if __name__ == "__main__":
    ut.main()
