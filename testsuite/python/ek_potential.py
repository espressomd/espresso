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
class EKPotential:
    BOX_L = [22.92, 28.65, 31.515]
    AGRID = 2.865
    TAU = 5.943635

    system = espressomd.System(box_l=BOX_L)
    system.time_step = TAU
    system.cell_system.skin = 0.4

    def tearDown(self):
        self.system.ekcontainer = None

    def test_potential(self):
        """
        Testing electrostatic potential of the EK
        """

        eps0 = 0.09135
        epsR = 8.0
        kT = 3.15379
        valency = 1.1
        density = 0.123

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

            ekspecies_pos[get_slice(0, dir)].density = density
            ekspecies_neg[get_slice(-1, dir)].density = density

            self.system.integrator.run(1)

            efield = - valency * density * self.AGRID**2 / \
                eps0 / epsR / self.BOX_L[dir]

            potential_simulation = np.copy(eksolver[:, :, :].potential)
            efield_simulation = np.gradient(
                potential_simulation, self.AGRID, axis=dir)

            atol = 1e-7 if self.lattice_params["single_precision"] else 2e-16

            np.testing.assert_allclose(
                efield_simulation, efield, rtol=0.0, atol=atol)

            in_plane_field = np.delete(np.arange(3), dir)
            for in_plane_dir in in_plane_field:
                efield_in_plane = np.gradient(
                    potential_simulation, self.AGRID, axis=in_plane_dir)
                np.testing.assert_allclose(
                    efield_in_plane, 0.0, rtol=0.0, atol=atol)


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class EKTestWalberla(EKPotential, ut.TestCase):
    lattice_params = {"single_precision": False, "gpu": False}


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class EKTestWalberlaSinglePrecision(EKPotential, ut.TestCase):
    lattice_params = {"single_precision": True, "gpu": False}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT", "CUDA"])
class EKTestWalberlaGPU(EKPotential, ut.TestCase):
    lattice_params = {"single_precision": False, "gpu": True}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT", "CUDA"])
class EKTestWalberlaSinglePrecisionGPU(EKPotential, ut.TestCase):
    lattice_params = {"single_precision": True, "gpu": True}


if __name__ == "__main__":
    ut.main()
