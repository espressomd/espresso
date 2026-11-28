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
class EKExternalElectricField:
    BOX_L = [17.19, 17.19, 17.19]
    AGRID = 2.865
    DIFFUSION_COEFFICIENT = 0.68
    TAU = 5.815
    TIMESTEPS = 10
    NUM_SAMPLES = 10

    np.random.seed(42)
    system = espressomd.System(box_l=BOX_L)
    system.time_step = TAU
    system.cell_system.skin = 0.4

    def tearDown(self):
        self.system.ekcontainer = None

    def test_external_electric_field(self):
        """
        Testing the flux caused by an external electric field.
        """

        kT = 3.64
        valency = 1.895
        external_electric_field = np.asarray([0.0, 0.0, 0.0])
        electric_field_max = 0.01

        density = 0.06

        lattice = espressomd.electrokinetics.Lattice(
            n_ghost_layers=1, agrid=self.AGRID)

        ekspecies_pos = espressomd.electrokinetics.EKSpecies(
            lattice=lattice, density=density, kT=kT, valency=valency,
            diffusion=self.DIFFUSION_COEFFICIENT, friction_coupling=False,
            advection=False, ext_efield=external_electric_field,
            tau=self.TAU, **self.lattice_params)
        ekspecies_neg = espressomd.electrokinetics.EKSpecies(
            lattice=lattice, density=density, kT=kT, valency=-valency,
            diffusion=self.DIFFUSION_COEFFICIENT, friction_coupling=False,
            advection=False, ext_efield=external_electric_field,
            tau=self.TAU, **self.lattice_params)

        eksolver = espressomd.electrokinetics.EKNone(
            lattice=lattice, tau=self.TAU, **self.lattice_params)

        self.system.ekcontainer = espressomd.electrokinetics.EKContainer(
            tau=self.TAU, solver=eksolver)
        self.system.ekcontainer.add(ekspecies_pos)
        self.system.ekcontainer.add(ekspecies_neg)

        atol = 5e-12 if not self.lattice_params["single_precision"] else 6e-10

        for _ in range(self.NUM_SAMPLES):
            external_electric_field = electric_field_max * np.random.random(3)
            ekspecies_pos.ext_efield = external_electric_field
            ekspecies_neg.ext_efield = external_electric_field
            self.system.integrator.run(self.TIMESTEPS)

            flux_pos = np.mean(np.copy(ekspecies_pos[:, :, :].flux), (0, 1, 2))
            flux_neg = np.mean(np.copy(ekspecies_neg[:, :, :].flux), (0, 1, 2))
            ref_flux = valency * density * \
                self.DIFFUSION_COEFFICIENT / kT * external_electric_field
            np.testing.assert_allclose(flux_pos, ref_flux, rtol=0.0, atol=atol)
            np.testing.assert_allclose(
                flux_neg, -ref_flux, rtol=0.0, atol=atol)


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class EKExternalElectricFieldDoublePrecisionCPU(EKExternalElectricField, ut.TestCase):
    lattice_params = {"single_precision": False, "gpu": False}


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class EKExternalElectricFieldSinglePrecisionCPU(EKExternalElectricField, ut.TestCase):
    lattice_params = {"single_precision": True, "gpu": False}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT", "CUDA"])
class EKExternalElectricFieldDoublePrecisionGPU(EKExternalElectricField, ut.TestCase):
    lattice_params = {"single_precision": False, "gpu": True}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT", "CUDA"])
class EKExternalElectricFieldSinglePrecisionGPU(EKExternalElectricField, ut.TestCase):
    lattice_params = {"single_precision": True, "gpu": True}


if __name__ == "__main__":
    ut.main()
