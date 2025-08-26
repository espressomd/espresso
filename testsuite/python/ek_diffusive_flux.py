#
# Copyright (C) 2022-2025 The ESPResSo project
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


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class EKDiffusiveFlux:
    BOX_L = [5., 5., 5.]
    AGRID = 1.0
    DIFFUSION_COEFFICIENT = 0.25
    TIMESTEPS = 1
    TAU = 1.0

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

        density = 1.0
        kT = 1.0

        lattice = self.ek_lattice_class(n_ghost_layers=1, agrid=self.AGRID)

        ekspecies = self.ek_species_class(
            lattice=lattice, density=0.0, kT=kT, valency=0.,
            diffusion=self.DIFFUSION_COEFFICIENT, friction_coupling=False,
            advection=False, tau=self.TAU, **self.ek_params)

        eksolver = espressomd.electrokinetics.EKNone(lattice=lattice)

        self.system.ekcontainer = espressomd.electrokinetics.EKContainer(
            tau=self.TAU, solver=eksolver)
        self.system.ekcontainer.add(ekspecies)

        center = np.array([2] * 3)
        ekspecies[center].density = density

        self.system.integrator.run(1)

        offset = np.array([-1, 0, 1])
        normalization_factor = 1.0 + 2. * np.sqrt(2) + 4.0 / 3.0 * np.sqrt(3)
        for x in offset:
            for y in offset:
                for z in offset:
                    direction = np.array([x, y, z])
                    local_flux = np.array(ekspecies[direction + center].flux)
                    dist = np.linalg.norm(direction)
                    if (dist > 0):
                        ref_flux = direction / normalization_factor * \
                            self.DIFFUSION_COEFFICIENT / dist / 2.0
                        np.testing.assert_allclose(
                            local_flux, ref_flux, rtol=1.0E-5)
                    else:
                        np.testing.assert_allclose(
                            local_flux, np.zeros(3), atol=1.0E-8)


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class EKTestWalberla(EKDiffusiveFlux, ut.TestCase):

    """Test for the waLBerla implementation of the EK in double-precision."""

    ek_lattice_class = espressomd.electrokinetics.LatticeWalberla
    ek_species_class = espressomd.electrokinetics.EKSpecies
    ek_params = {"single_precision": False}


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class EKTestWalberlaSinglePrecision(EKDiffusiveFlux, ut.TestCase):

    """Test for the waLBerla implementation of the EK in single-precision."""

    ek_lattice_class = espressomd.electrokinetics.LatticeWalberla
    ek_species_class = espressomd.electrokinetics.EKSpecies
    ek_params = {"single_precision": True}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT", "CUDA"])
class EKTestWalberlaGPU(EKDiffusiveFlux, ut.TestCase):

    """Test for the waLBerla implementation of the EK in double-precision GPU."""

    ek_lattice_class = espressomd.electrokinetics.LatticeWalberla
    ek_species_class = espressomd.electrokinetics.EKSpeciesGPU
    ek_params = {"single_precision": False}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT", "CUDA"])
class EKTestWalberlaSinglePrecisionGPU(EKDiffusiveFlux, ut.TestCase):

    """Test for the waLBerla implementation of the EK in single-precision GPU."""

    ek_lattice_class = espressomd.electrokinetics.LatticeWalberla
    ek_species_class = espressomd.electrokinetics.EKSpeciesGPU
    ek_params = {"single_precision": True}


if __name__ == "__main__":
    ut.main()
