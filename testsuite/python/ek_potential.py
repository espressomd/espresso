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
class EKEOF:
    BOX_L = [15., 18., 21.]
    AGRID = 1.5
    DIFFUSION_COEFFICIENT = 0.0
    TAU = 1.6

    system = espressomd.System(box_l=BOX_L)
    system.time_step = TAU
    system.cell_system.skin = 0.4

    def tearDown(self):
        self.system.ekcontainer = None

    def test_potential(self):
        """
        Testing electrostatic potential of the EK
        """

        eps0 = 0.015
        epsR = 18.5
        kT = 2.
        valency = 1.1

        density = 0.0
        external_electric_field = np.asarray([0.0, 0.0, 0.0])

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

        # Simple linear potential test
        def get_slice(value, dir):
            match dir:
                case 0:
                    return (value, slice(None), slice(None))
                case 1:
                    return (slice(None), value, slice(None))
                case 2:
                    return (slice(None), slice(None), value)

        for dir in range(3):
            ekspecies_pos[:, :, :].density = 0.0
            ekspecies_neg[:, :, :].density = 0.0

            ekspecies_pos[get_slice(0, dir)].density = 1.0
            ekspecies_neg[get_slice(-1, dir)].density = 1.0

            self.system.integrator.run(1)

            pot_min = eksolver[0, 0, 0].potential  # fixing the offset
            self.assertTrue(np.isfinite(pot_min), msg="Potential grid is NaN")

            # We divide by BOX_L, because over the PBC the surfaces are
            # one unit apart, which defines the voltage
            ref_voltage = self.AGRID * valency * \
                (self.BOX_L[dir] - self.AGRID) / eps0 / epsR / self.BOX_L[dir]
            width = int(self.BOX_L[dir] / self.AGRID)
            slope = -ref_voltage / (width - 1)

            for d in range(width):
                ref_value = d * slope + pot_min
                potential = np.copy(
                    eksolver[get_slice(d, dir)].potential).flatten()
                np.testing.assert_allclose(
                    potential, ref_value, rtol=self.rtol)


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class EKTestWalberla(EKEOF, ut.TestCase):

    """Test for the waLBerla implementation of the EK in double-precision."""

    ek_lattice_class = espressomd.electrokinetics.Lattice
    ek_species_class = espressomd.electrokinetics.EKSpecies
    ek_solver_class = espressomd.electrokinetics.EKFFT
    ek_params = {"single_precision": False, "gpu": False}
    rtol = 5e-12


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class EKTestWalberlaSinglePrecision(EKEOF, ut.TestCase):

    """Test for the waLBerla implementation of the EK in single-precision."""

    ek_lattice_class = espressomd.electrokinetics.Lattice
    ek_species_class = espressomd.electrokinetics.EKSpecies
    ek_solver_class = espressomd.electrokinetics.EKFFT
    ek_params = {"single_precision": True, "gpu": False}
    rtol = 1e-5


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT", "CUDA"])
class EKTestWalberlaGPU(EKEOF, ut.TestCase):

    """Test for the waLBerla implementation of the EK in double-precision."""

    ek_lattice_class = espressomd.electrokinetics.Lattice
    ek_species_class = espressomd.electrokinetics.EKSpecies
    ek_solver_class = espressomd.electrokinetics.EKFFT
    ek_params = {"single_precision": False, "gpu": True}
    rtol = 5e-12


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT", "CUDA"])
class EKTestWalberlaSinglePrecisionGPU(EKEOF, ut.TestCase):

    """Test for the waLBerla implementation of the EK in single-precision."""

    ek_lattice_class = espressomd.electrokinetics.Lattice
    ek_species_class = espressomd.electrokinetics.EKSpecies
    ek_solver_class = espressomd.electrokinetics.EKFFT
    ek_params = {"single_precision": True, "gpu": True}
    rtol = 1e-5


if __name__ == "__main__":
    ut.main()
