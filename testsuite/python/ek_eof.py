#
# Copyright (C) 2022-2023 The ESPResSo project
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
import scipy.optimize
import unittest as ut
import unittest_decorators as utx

import espressomd
import espressomd.lb
import espressomd.shapes
import espressomd.electrokinetics


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class EKEOF:
    BOX_L = [45., 9., 9.]
    AGRID = 1.5
    DENSITY = 1
    DIFFUSION_COEFFICIENT = 0.25
    TIMESTEPS = 5000
    TAU = 1.6

    system = espressomd.System(box_l=BOX_L)
    system.time_step = TAU
    system.cell_system.skin = 0.4

    def tearDown(self):
        self.system.actors.clear()
        self.system.ekcontainer.clear()

    def test_eof(self):
        """
        Testing EK for the electroosmotic flow
        """

        eps0 = 0.015
        epsR = 18.5
        kT = 2.
        offset = self.AGRID
        d = self.system.box_l[0] - 2 * offset
        valency = 1.1
        external_electric_field = np.asarray([0.0, 0.001, 0.0])

        visc = 1. / 6.
        eta = 1.0 * visc

        density = 0.0006

        lattice = self.ek_lattice_class(n_ghost_layers=1, agrid=self.AGRID)

        ekspecies = self.ek_species_class(
            lattice=lattice, density=density, kT=kT, valency=valency,
            diffusion=self.DIFFUSION_COEFFICIENT, friction_coupling=True,
            advection=True, ext_efield=external_electric_field,
            tau=self.TAU, **self.ek_params)
        ekwallcharge = self.ek_species_class(
            lattice=lattice, density=0.0, kT=kT, diffusion=0.0, tau=self.TAU,
            valency=-valency, friction_coupling=False, advection=False,
            **self.ek_params)

        eksolver = self.ek_solver_class(
            lattice=lattice, permittivity=eps0 * epsR, **self.ek_params)

        self.system.ekcontainer.tau = self.TAU
        self.system.ekcontainer.solver = eksolver
        self.system.ekcontainer.add(ekspecies)
        self.system.ekcontainer.add(ekwallcharge)

        lb_fluid = espressomd.lb.LBFluidWalberla(
            lattice=lattice, density=1.0, kinematic_viscosity=visc,
            tau=self.TAU, **self.ek_params)
        self.system.actors.add(lb_fluid)

        wall_bot = espressomd.shapes.Wall(normal=[1, 0, 0], dist=offset)
        wall_top = espressomd.shapes.Wall(
            normal=[-1, 0, 0], dist=-self.BOX_L[0] + offset)
        for obj in (wall_bot, wall_top):
            ekspecies.add_boundary_from_shape(
                obj, [0.0, 0.0, 0.0], espressomd.electrokinetics.FluxBoundary)
            ekspecies.add_boundary_from_shape(
                obj, 0.0, espressomd.electrokinetics.DensityBoundary)
            lb_fluid.add_boundary_from_shape(obj, [0.0, 0.0, 0.0])

        ekspecies[0, :, :].density = 0.0
        ekspecies[-1, :, :].density = 0.0

        density_wall = density * d / 2
        sigma = -valency * density_wall

        ekwallcharge[0, :, :].density = density_wall
        ekwallcharge[-1, :, :].density = density_wall

        self.system.integrator.run(self.TIMESTEPS)

        def transcendental_equation(
                x, valency, distance, kT, sigma, eps0, epsR):
            return x * np.tan(valency * distance / (4. * kT) * x) + \
                sigma / (eps0 * epsR)

        integration_constant = scipy.optimize.fsolve(
            lambda x: transcendental_equation(
                x=x,
                valency=valency,
                distance=d,
                kT=kT,
                sigma=sigma,
                eps0=eps0,
                epsR=epsR),
            0.1)

        def calc_analytic_density(
                x, integration_constant, valency, kT, eps0, epsR):
            return (epsR * eps0) * integration_constant**2 / (2 * kT) / np.square(
                np.cos(valency * integration_constant / (2 * kT) * x))

        def calc_analytic_velocity(x, integration_constant, valency,
                                   distance, kT, eps0, epsR, eta, external_electric_field):
            return 2 * np.linalg.norm(external_electric_field) * epsR * eps0 * kT / (eta * valency) * (
                np.log(np.cos(valency * integration_constant / (2 * kT) * x)) - np.log(np.cos(distance * valency * integration_constant / (4 * kT))))

        simulated_density = ekspecies[:, 0, 0].density.squeeze()
        simulated_velocity = lb_fluid[:, 0, 0].velocity.squeeze()[:, 1]

        x_sim = (np.arange(lattice.shape[0]) -
                 lattice.shape[0] / 2 + 0.5) * self.AGRID

        analytic_density = calc_analytic_density(
            x=x_sim,
            integration_constant=integration_constant,
            valency=valency,
            kT=kT,
            eps0=eps0,
            epsR=epsR)
        analytic_density[np.logical_or(
            x_sim < -self.system.box_l[0] / 2 + offset,
            x_sim > self.system.box_l[0] / 2 - offset)] = 0.
        np.testing.assert_allclose(
            simulated_density, analytic_density, rtol=2e-2)

        analytic_velocity = calc_analytic_velocity(
            x=x_sim,
            integration_constant=integration_constant,
            valency=valency,
            distance=d,
            kT=kT,
            eps0=eps0,
            epsR=epsR,
            eta=eta,
            external_electric_field=external_electric_field)
        analytic_velocity[np.logical_or(
            x_sim < -self.system.box_l[0] / 2 + offset,
            x_sim > self.system.box_l[0] / 2 - offset)] = 0.
        np.testing.assert_allclose(
            simulated_velocity,
            analytic_velocity,
            rtol=2e-2)


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class EKTestWalberla(EKEOF, ut.TestCase):

    """Test for the Walberla implementation of the EK in double-precision."""

    ek_lattice_class = espressomd.electrokinetics.LatticeWalberla
    ek_species_class = espressomd.electrokinetics.EKSpecies
    ek_solver_class = espressomd.electrokinetics.EKFFT
    ek_params = {"single_precision": False}


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class EKTestWalberlaSinglePrecision(EKEOF, ut.TestCase):

    """Test for the Walberla implementation of the EK in single-precision."""

    ek_lattice_class = espressomd.electrokinetics.LatticeWalberla
    ek_species_class = espressomd.electrokinetics.EKSpecies
    ek_solver_class = espressomd.electrokinetics.EKFFT
    ek_params = {"single_precision": True}


if __name__ == "__main__":
    ut.main()
