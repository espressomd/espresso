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


@utx.skipIfMissingFeatures(["WALBERLA"])
class EKFixedDensity(ut.TestCase):
    BOX_L = 42.
    AGRID = 1.0
    DENSITY = 1
    DIFFUSION_COEFFICIENT = 0.5
    TIME = 5000

    INLET_CONCENTRATION = 1.0
    OUTLET_CONCENTRATION = 0.01

    system = espressomd.System(box_l=[BOX_L, 3, 3])
    system.time_step = 1.0
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

        lattice = espressomd.lb.LatticeWalberla(
            n_ghost_layers=1, agrid=self.AGRID)

        ekspecies = espressomd.EKSpecies.EKSpecies(
            lattice=lattice, density=0.0, diffusion=self.DIFFUSION_COEFFICIENT,
            kT=0.0, valency=0.0, advection=False, friction_coupling=False,
            ext_efield=[0, 0, 0], single_precision=single_precision)

        eksolver = espressomd.EKSpecies.EKNone(lattice=lattice)

        self.system.ekcontainer.add(ekspecies)
        self.system.ekcontainer.tau = 1.0
        self.system.ekcontainer.solver = eksolver

        # left and right no flux
        ekspecies[0, :, :].flux_boundary = \
            espressomd.EKSpecies.FluxBoundary([0, 0, 0])
        ekspecies[-1, :, :].flux_boundary = \
            espressomd.EKSpecies.FluxBoundary([0, 0, 0])

        left_slice = ekspecies[1, :, :]
        left_slice.density = 1.0
        left_slice.density_boundary = espressomd.EKSpecies.DensityBoundary(
            self.INLET_CONCENTRATION)

        right_slice = ekspecies[-2, :, :]
        right_slice.density_boundary = espressomd.EKSpecies.DensityBoundary(
            self.OUTLET_CONCENTRATION)

        self.system.integrator.run(self.TIME)

        effective_boxl = self.BOX_L - 2
        domain_positions = np.arange(effective_boxl, dtype=np.float64)

        measured_values = ekspecies[1:-1, 1, 1].density.squeeze()

        slope = (self.OUTLET_CONCENTRATION -
                 self.INLET_CONCENTRATION) / (effective_boxl - 1)
        offset = self.INLET_CONCENTRATION
        analytic_values = slope * domain_positions + offset

        np.testing.assert_almost_equal(
            measured_values,
            analytic_values,
            decimal_precision)


if __name__ == "__main__":
    ut.main()
