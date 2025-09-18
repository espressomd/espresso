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
class EKIonConductivity:
    BOX_L = [9., 9., 9.]
    AGRID = 1.5
    DIFFUSION_COEFFICIENT = 0.25
    TIMESTEPS = 10
    TAU = 1.6
    NUM_SAMPLES = 10

    np.random.seed(37)
    system = espressomd.System(box_l=BOX_L)
    system.time_step = TAU
    system.cell_system.skin = 0.4

    def tearDown(self):
        self.system.ekcontainer = None

    def test_conductivity(self):
        """
        Testing the ion conductivity of a ionic solution by measuring
        the flux.
        """

        eps0 = 0.015
        epsR = 18.5
        kT = 2.
        valency = 1.1
        external_electric_field = np.asarray([0.0, 0.0, 0.0])
        electric_field_max = 0.01

        density = 0.0006

        lattice = self.ek_lattice_class(n_ghost_layers=1, agrid=self.AGRID)

        ekspecies_pos = self.ek_species_class(
            lattice=lattice, density=density, kT=kT, valency=valency,
            diffusion=self.DIFFUSION_COEFFICIENT, friction_coupling=False,
            advection=False, ext_efield=external_electric_field,
            tau=self.TAU, **self.ek_params)
        ekspecies_neg = self.ek_species_class(
            lattice=lattice, density=density, kT=kT, valency=-valency,
            diffusion=self.DIFFUSION_COEFFICIENT, friction_coupling=False,
            advection=False, ext_efield=external_electric_field,
            tau=self.TAU, **self.ek_params)

        eksolver = self.ek_solver_class(
            lattice=lattice, permittivity=eps0 * epsR, **self.ek_params)

        self.system.ekcontainer = espressomd.electrokinetics.EKContainer(
            tau=self.TAU, solver=eksolver)
        self.system.ekcontainer.add(ekspecies_pos)
        self.system.ekcontainer.add(ekspecies_neg)

        for _ in range(self.NUM_SAMPLES):
            external_electric_field = electric_field_max * np.random.random(3)
            ekspecies_pos.ext_efield = external_electric_field
            ekspecies_neg.ext_efield = external_electric_field
            self.system.integrator.run(self.TIMESTEPS)

            conductivity = np.mean(
                ekspecies_pos[:, :, :].flux - ekspecies_neg[:, :, :].flux, (0, 1, 2))
            conductivity /= external_electric_field
            ref_value = 2 * valency * density * self.DIFFUSION_COEFFICIENT / kT
            np.testing.assert_allclose(conductivity, ref_value, rtol=5e-4)


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class EKTestWalberla(EKIonConductivity, ut.TestCase):

    """Test for the waLBerla implementation of the EK in double-precision."""

    ek_lattice_class = espressomd.electrokinetics.LatticeWalberla
    ek_species_class = espressomd.electrokinetics.EKSpecies
    ek_solver_class = espressomd.electrokinetics.EKFFT
    ek_params = {"single_precision": False}


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class EKTestWalberlaSinglePrecision(EKIonConductivity, ut.TestCase):

    """Test for the waLBerla implementation of the EK in single-precision."""

    ek_lattice_class = espressomd.electrokinetics.LatticeWalberla
    ek_species_class = espressomd.electrokinetics.EKSpecies
    ek_solver_class = espressomd.electrokinetics.EKFFT
    ek_params = {"single_precision": True}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT", "CUDA"])
class EKTestWalberlaGPU(EKIonConductivity, ut.TestCase):

    """Test for the waLBerla implementation of the EK in double-precision GPU."""

    ek_lattice_class = espressomd.electrokinetics.LatticeWalberla
    ek_species_class = espressomd.electrokinetics.EKSpeciesGPU
    ek_solver_class = espressomd.electrokinetics.EKFFTGPU
    ek_params = {"single_precision": False}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT", "CUDA"])
class EKTestWalberlaSinglePrecisionGPU(EKIonConductivity, ut.TestCase):

    """Test for the waLBerla implementation of the EK in single-precision GPU."""

    ek_lattice_class = espressomd.electrokinetics.LatticeWalberla
    ek_species_class = espressomd.electrokinetics.EKSpeciesGPU
    ek_solver_class = espressomd.electrokinetics.EKFFTGPU
    ek_params = {"single_precision": True}


if __name__ == "__main__":
    ut.main()
