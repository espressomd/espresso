#
# Copyright (C) 2022 The ESPResSo project
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
import numpy as np
import espressomd.shapes
import scipy.optimize


@utx.skipIfMissingFeatures(["WALBERLA", "WALBERLA_FFT"])
class EKEOF(ut.TestCase):
    BOX_L = [30., 6., 6.]
    AGRID = 1.0
    DENSITY = 1
    DIFFUSION_COEFFICIENT = 0.25
    TIME = 3000

    system = espressomd.System(box_l=BOX_L)
    system.time_step = 1.0
    system.cell_system.skin = 0.4

    def tearDown(self) -> None:
        self.system.actors.clear()
        self.system.ekcontainer.clear()

    def analytical_density(pos: np.ndarray, time: int, D: float):
        return (4 * np.pi * D * time)**(-3 / 2) * \
            np.exp(-np.sum(np.square(pos), axis=-1) / (4 * D * time))

    def test_eof_single(self):
        self.detail_test_eof(single_precision=True)

    def test_eof_double(self):
        self.detail_test_eof(single_precision=False)

    def detail_test_eof(self, single_precision: bool):
        """
        Testing EK for the electroosmotic flow
        """

        eps0 = 1.0
        epsR = 1.0
        kT = 1.0
        offset = 1
        d = self.system.box_l[0] - 2 * offset
        valency = 1.0
        external_electric_field = np.asarray([0.0, 0.001, 0.0])

        visc = 1. / 6.
        eta = 1.0 * visc

        density = 0.0006

        lattice = espressomd.lb.LatticeWalberla(
            n_ghost_layers=1, agrid=self.AGRID)

        ekspecies = espressomd.EKSpecies.EKSpecies(lattice=lattice,
                                                   density=density, kT=kT, diffusion=self.DIFFUSION_COEFFICIENT, valency=valency,
                                                   advection=True, friction_coupling=True, ext_efield=external_electric_field, single_precision=single_precision)
        ekwallcharge = espressomd.EKSpecies.EKSpecies(lattice=lattice,
                                                      density=0.0, kT=kT, diffusion=0.0, valency=-valency, advection=False, friction_coupling=False, ext_efield=[0, 0, 0], single_precision=single_precision)

        eksolver = espressomd.EKSpecies.EKFFT(
            lattice=lattice, permittivity=eps0 * epsR, single_precision=single_precision)
        self.system.ekcontainer.add(ekspecies)
        self.system.ekcontainer.add(ekwallcharge)

        self.system.ekcontainer.tau = 1.0
        self.system.ekcontainer.solver = eksolver

        lb_fluid = espressomd.lb.LBFluidWalberla(
            lattice=lattice, density=1.0, viscosity=visc, tau=1.0, single_precision=single_precision)
        self.system.actors.add(lb_fluid)

        wall_bot = espressomd.shapes.Wall(normal=[1, 0, 0], dist=offset)
        wall_top = espressomd.shapes.Wall(
            normal=[-1, 0, 0], dist=-self.BOX_L[0] + offset)
        for obj in (wall_bot, wall_top):
            ekspecies.add_boundary_from_shape(
                obj, [0.0, 0.0, 0.0], espressomd.EKSpecies.FluxBoundary)
            lb_fluid.add_boundary_from_shape(obj, [0.0, 0.0, 0.0])

        ekspecies[0, :, :].density = 0.0
        ekspecies[-1, :, :].density = 0.0

        sigma = density * d / 2

        ekwallcharge[0, :, :].density = sigma
        ekwallcharge[-1, :, :].density = sigma

        self.system.integrator.run(self.TIME)

        # TODO: figure out if the runtime is enough
        # TODO: make an actual comparison with the analytic solution
        # TODO: add dielectric permittivities? electrostatic prefactor? to
        # EKFFT

        def trans(x):
            return x * np.tan(valency * d / (4. * kT) * x) - \
                sigma / (eps0 * epsR)

        C = scipy.optimize.fsolve(trans, 0.1)

        def analytical(x):
            return (epsR * eps0) * C ** 2 / (2 * kT) / np.square(
                np.cos(valency * C / (2 * kT) * x))

        def vel(x):
            return 2 * np.linalg.norm(external_electric_field) * epsR * eps0 * kT / (eta * valency) * (
                np.log(np.cos(valency * C / (2 * kT) * x)) - np.log(np.cos(d * valency * C / (4 * kT))))

        # import matplotlib.pyplot as plt
        #
        # plt.figure()

        # simulated_density = ekspecies[:, 0, 0].density.squeeze()
        simulated_velocity = lb_fluid[:, 0, 0].velocity.squeeze()[:, 1]

        x_sim = np.arange(
            self.system.box_l[0]) - self.system.box_l[0] / 2 + 0.5
        # x_shifted = np.linspace(x_sim[0], x_sim[-1], 200)
        # ions boundary is half a cell offset to the fluid cell
        # x_ions = np.linspace(x_sim[0] - 0.5, x_sim[-1] + 0.5, 200)
        # plt.plot(x_sim, simulated_density, "x", label="density")
        # plt.plot(x_sim, simulated_velocity, "x", label="y-velocity")

        analytic_density = analytical(x_sim)
        analytic_density[np.logical_or(
            x_sim < -self.system.box_l[0] / 2 + offset,
            x_sim > self.system.box_l[0] / 2 - offset)] = 0.
        # np.testing.assert_allclose(
        #     simulated_density, analytic_density, rtol=2e-3)

        # TODO: figure out how to test the velocity profile...z
        analytic_velocity = vel(x_sim)
        analytic_velocity[np.logical_or(
            x_sim < -self.system.box_l[0] / 2 + offset,
            x_sim > self.system.box_l[0] / 2 - offset)] = 0.
        np.testing.assert_allclose(
            simulated_velocity,
            analytic_velocity,
            rtol=3e-2)
        # print((simulated_velocity - analytic_velocity) / analytic_velocity)
        # print(np.sum(np.abs(simulated_velocity - analytic_velocity)))

        # analytic_density = analytical(x_ions)
        # analytic_velocity = vel(x_shifted)
        # analytic_density[np.logical_or(x_ions < - self.system.box_l[0] / 2 + offset - 0.5, x_ions > self.system.box_l[0] / 2 - offset + 0.5)] = 0.
        # analytic_velocity[np.logical_or(x_shifted < - self.system.box_l[0] / 2 + offset, x_shifted > self.system.box_l[0] / 2 - offset)] = 0.

        # plt.plot(x_ions, analytic_density, label="analytic density")
        # plt.plot(x_shifted, analytic_velocity, label="analytic velocity")

        # plt.legend()

        # plt.show()


if __name__ == "__main__":
    ut.main()
