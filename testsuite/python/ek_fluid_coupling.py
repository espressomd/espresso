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
import espressomd.lb
import espressomd.shapes
import espressomd.electrokinetics


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class Fluid_coupling:
    BOX_L = [9., 9., 9.]
    AGRID = 1.5
    DENSITY = 1
    DIFFUSION_COEFFICIENT = 0.25
    TIMESTEPS = 50
    TAU = 1.6

    system = espressomd.System(box_l=BOX_L)
    system.time_step = TAU
    system.cell_system.skin = 0.4

    def tearDown(self):
        self.system.lb = None
        self.system.ekcontainer = None

    def test_fluid_coupling(self):
        """
        Testing the friction forces acting on the LB for multiple EK species
        """
        for opposite_external_fields in [True, False]:
            with self.subTest(msg=f"{opposite_external_fields=}"):
                self.check_fluid_coupling(opposite_external_fields)
                self.tearDown()

    def check_fluid_coupling(self, opposite_external_fields):
        eps0 = 0.015
        epsR = 18.5
        kT = 2.
        fluid_density = 1.2
        valency = 1.1
        external_electric_field_pos = (
            np.random.rand(3) - np.array([0.5] * 3)) * 0.01

        if opposite_external_fields:
            external_electric_field_neg = -external_electric_field_pos
        else:
            external_electric_field_neg = external_electric_field_pos

        visc = 1. / 6.

        density = 0.0006

        lattice = self.ek_lattice_class(n_ghost_layers=2, agrid=self.AGRID)

        ekspecies_pos = self.ek_species_class(
            lattice=lattice, density=density, kT=kT, valency=valency,
            diffusion=self.DIFFUSION_COEFFICIENT, friction_coupling=True,
            advection=True, ext_efield=external_electric_field_pos,
            tau=self.TAU, **self.ek_params)

        ekspecies_neg = self.ek_species_class(
            lattice=lattice, density=density, kT=kT, valency=-valency,
            diffusion=self.DIFFUSION_COEFFICIENT, friction_coupling=True,
            advection=True, ext_efield=external_electric_field_neg,
            tau=self.TAU, **self.ek_params)

        eksolver = self.ek_solver_class(
            lattice=lattice, permittivity=eps0 * epsR, **self.ek_params)

        self.system.ekcontainer = espressomd.electrokinetics.EKContainer(
            tau=self.TAU, solver=eksolver)
        self.system.ekcontainer.add(ekspecies_pos)
        self.system.ekcontainer.add(ekspecies_neg)

        lb_fluid = self.lb_class(
            lattice=lattice, density=fluid_density, kinematic_viscosity=visc,
            tau=self.TAU, **self.lb_params)
        self.system.lb = lb_fluid

        self.system.integrator.run(self.TIMESTEPS)

        forces = np.copy(lb_fluid[:, :, :].last_applied_force)

        if opposite_external_fields:
            expected_force = 3 * external_electric_field_pos * \
                valency * density * self.AGRID**2
            expected_force = np.full_like(forces, expected_force)
            rtol = 1e-4 if self.ek_params["single_precision"] else 1e-5
            np.testing.assert_allclose(
                forces, expected_force, rtol=rtol, atol=1e-10)
        else:
            np.testing.assert_allclose(
                forces, np.zeros_like(forces), atol=1e-10)


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class EKTestWalberlaDoublePrecisionCPU(Fluid_coupling, ut.TestCase):

    """Test for the Walberla implementation of the EK in double-precision."""

    ek_lattice_class = espressomd.electrokinetics.Lattice
    ek_species_class = espressomd.electrokinetics.EKSpecies
    ek_solver_class = espressomd.electrokinetics.EKFFT
    lb_class = espressomd.lb.LBFluid
    lb_params = {"single_precision": False, "gpu": False}
    ek_params = {"single_precision": False, "gpu": False}


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class EKTestWalberlaSinglePrecisionCPU(Fluid_coupling, ut.TestCase):

    """Test for the Walberla implementation of the EK in single-precision."""

    ek_lattice_class = espressomd.electrokinetics.Lattice
    ek_species_class = espressomd.electrokinetics.EKSpecies
    ek_solver_class = espressomd.electrokinetics.EKFFT
    lb_class = espressomd.lb.LBFluid
    lb_params = {"single_precision": True, "gpu": False}
    ek_params = {"single_precision": True, "gpu": False}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT", "CUDA"])
class EKTestWalberlaDoublePrecisionGPU(Fluid_coupling, ut.TestCase):

    """Test for the Walberla implementation of the EK in double-precision."""

    ek_lattice_class = espressomd.electrokinetics.Lattice
    ek_species_class = espressomd.electrokinetics.EKSpecies
    ek_solver_class = espressomd.electrokinetics.EKFFT
    lb_class = espressomd.lb.LBFluid
    lb_params = {"single_precision": False, "gpu": True}
    ek_params = {"single_precision": False, "gpu": True}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT", "CUDA"])
class EKTestWalberlaSinglePrecisionGPU(Fluid_coupling, ut.TestCase):

    """Test for the Walberla implementation of the EK in single-precision."""

    ek_lattice_class = espressomd.electrokinetics.Lattice
    ek_species_class = espressomd.electrokinetics.EKSpecies
    ek_solver_class = espressomd.electrokinetics.EKFFT
    lb_class = espressomd.lb.LBFluid
    lb_params = {"single_precision": True, "gpu": True}
    ek_params = {"single_precision": True, "gpu": True}


if __name__ == "__main__":
    ut.main()
