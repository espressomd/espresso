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
import unittest as ut
import unittest_decorators as utx

import espressomd
import espressomd.electrokinetics


@utx.skipIfMissingFeatures(["WALBERLA"])
class EKFixedDensity(ut.TestCase):
    AGRID = 1.1
    BOX_L = np.array([32., 3., 3.]) * AGRID
    DENSITY = 1

    DIFFUSION_COEFFICIENT = 0.2
    TIME = 5000
    TAU = 1.81

    INLET_CONCENTRATION = 1.0
    OUTLET_CONCENTRATION = 0.01

    system = espressomd.System(box_l=BOX_L)
    system.time_step = TAU
    system.cell_system.skin = 0.4

    def tearDown(self) -> None:
        self.system.ekcontainer.clear()

    def test_constant_density_bc_single(self):
        self.detail_test_constant_density_bc(single_precision=True)

    def test_constant_density_bc_double(self):
        self.detail_test_constant_density_bc(single_precision=False)

    def detail_test_constant_density_bc(self, single_precision: bool):
        """ effective 1D system with linear equilibrium profile """

        decimal_precision: int = 5 if single_precision else 7

        lattice = espressomd.electrokinetics.LatticeWalberla(
            n_ghost_layers=1, agrid=self.AGRID)

        ekspecies = espressomd.electrokinetics.EKSpecies(
            lattice=lattice, density=0.0, diffusion=self.DIFFUSION_COEFFICIENT,
            valency=0.0, advection=False, friction_coupling=False,
            single_precision=single_precision, tau=self.TAU)

        eksolver = espressomd.electrokinetics.EKNone(lattice=lattice)

        self.system.ekcontainer.tau = self.TAU
        self.system.ekcontainer.solver = eksolver
        self.system.ekcontainer.add(ekspecies)

        # left and right no flux
        ekspecies[0, :, :].flux_boundary = \
            espressomd.electrokinetics.FluxBoundary([0, 0, 0])
        ekspecies[-1, :, :].flux_boundary = \
            espressomd.electrokinetics.FluxBoundary([0, 0, 0])

        left_slice = ekspecies[1, :, :]
        left_slice.density = 1.0
        left_slice.density_boundary = espressomd.electrokinetics.DensityBoundary(
            self.INLET_CONCENTRATION)

        right_slice = ekspecies[-2, :, :]
        right_slice.density_boundary = espressomd.electrokinetics.DensityBoundary(
            self.OUTLET_CONCENTRATION)

        self.system.integrator.run(self.TIME)

        effective_boxl = (lattice.shape[0] - 3) * self.AGRID
        domain_positions = np.arange(
            lattice.shape[0] - 2,
            dtype=np.float64) * self.AGRID

        measured_values = ekspecies[1:-1, 1, 1].density.squeeze()

        slope = (self.OUTLET_CONCENTRATION -
                 self.INLET_CONCENTRATION) / effective_boxl
        offset = self.INLET_CONCENTRATION
        analytic_values = slope * domain_positions + offset

        np.testing.assert_almost_equal(
            measured_values,
            analytic_values,
            decimal_precision)


if __name__ == "__main__":
    ut.main()
