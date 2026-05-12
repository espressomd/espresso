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
import espressomd.electrokinetics


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class EKDiffusiveFlux:
    BOX_L = [14.325, 14.325, 14.325]
    AGRID = 2.865
    DENSITY = 1
    DIFFUSION_COEFFICIENT = 0.68
    TIMESTEPS = 1
    TAU = 1.815

    system = espressomd.System(box_l=BOX_L)
    system.time_step = TAU
    system.cell_system.skin = 0.4

    def tearDown(self):
        self.system.ekcontainer = None

    def test_diffusive_flux(self):
        """
        Testing the ion conductivity of a ionic solution by measuring
        the flux.
        """

        kT = 19.332

        lattice = espressomd.electrokinetics.Lattice(
            n_ghost_layers=1, agrid=self.AGRID)

        ekspecies = espressomd.electrokinetics.EKSpecies(
            lattice=lattice, density=0.0, kT=kT, valency=0.,
            diffusion=self.DIFFUSION_COEFFICIENT, friction_coupling=False,
            advection=False, tau=self.TAU, **self.lattice_params)

        eksolver = espressomd.electrokinetics.EKNone(
            lattice=lattice, tau=self.TAU)

        self.system.ekcontainer = espressomd.electrokinetics.EKContainer(
            tau=self.TAU, solver=eksolver)
        self.system.ekcontainer.add(ekspecies)

        center = np.array([2] * 3)
        ekspecies[center].density = 1.0

        self.system.integrator.run(1)

        atol = 3e-9 if self.lattice_params["single_precision"] else 7e-12

        offset = np.array([-1, 0, 1])
        normalization_factor = 1.0 + 2. * np.sqrt(2) + 4.0 / 3.0 * np.sqrt(3)
        for x in offset:
            for y in offset:
                for z in offset:
                    direction = np.array([x, y, z])
                    local_flux = np.array(ekspecies[direction + center].flux)
                    dist = np.linalg.norm(direction)
                    if dist > 0:
                        ref_flux = direction / normalization_factor * \
                            self.DIFFUSION_COEFFICIENT / dist / 2.0 / self.AGRID
                        np.testing.assert_allclose(
                            local_flux, ref_flux, rtol=0.0, atol=atol)
                    else:
                        np.testing.assert_allclose(
                            local_flux, np.zeros(3), rtol=0.0, atol=2 * atol)


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class EKTestWalberlaDoublePrecisionCPU(EKDiffusiveFlux, ut.TestCase):
    lattice_params = {"single_precision": False, "gpu": False}


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class EKTestWalberlaSinglePrecisionCPU(EKDiffusiveFlux, ut.TestCase):
    lattice_params = {"single_precision": True, "gpu": False}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT", "CUDA"])
class EKTestWalberlaDoublePrecisionGPU(EKDiffusiveFlux, ut.TestCase):
    lattice_params = {"single_precision": False, "gpu": True}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT", "CUDA"])
class EKTestWalberlaSinglePrecisionGPU(EKDiffusiveFlux, ut.TestCase):
    lattice_params = {"single_precision": True, "gpu": True}


if __name__ == "__main__":
    ut.main()
