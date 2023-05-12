#
# Copyright (C) 2010-2022 The ESPResSo project
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

import espressomd
import espressomd.lb
import unittest as ut
import unittest_decorators as utx
import numpy as np

# Define the LB parameters
TIME_STEP = 0.008
AGRID = .4
GRID_SIZE = 6
KVISC = 4
DENS = 2.3
F = 5.5 / GRID_SIZE**3
GAMMA = 1


LB_PARAMS = {'agrid': AGRID,
             'density': DENS,
             'kinematic_viscosity': KVISC,
             'tau': TIME_STEP,
             'ext_force_density': np.array([-.7 * F, .9 * F, .8 * F])}


class TestLBMomentumConservation:
    """
    Tests momentum conservation for an LB coupled to a particle, where opposing
    forces are applied to LB and particle. The test should uncover issues
    with boundary and ghost layer handling.

    """

    system = espressomd.System(box_l=[GRID_SIZE * AGRID] * 3)
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.01
    n_nodes = system.cell_system.get_state()["n_nodes"]

    def setUp(self):
        self.set_cellsystem()
        self.lbf = self.lb_class(**LB_PARAMS, **self.lb_params)

    def tearDown(self):
        self.system.actors.clear()
        self.system.thermostat.turn_off()
        self.system.part.clear()

    def test(self):
        self.system.actors.add(self.lbf)
        self.system.thermostat.set_lb(LB_fluid=self.lbf, gamma=GAMMA, seed=1)
        np.testing.assert_allclose(
            self.lbf.ext_force_density,
            LB_PARAMS["ext_force_density"])

        # Initial momentum before integration = 0
        mom_tol = 1E-4 if self.lbf.single_precision else 1E-12
        np.testing.assert_allclose(
            self.system.analysis.linear_momentum(), [0., 0., 0.], atol=mom_tol)

        ext_fluid_force = self.system.volume() * LB_PARAMS["ext_force_density"]

        p = self.system.part.add(
            pos=self.system.box_l / 2, ext_force=-ext_fluid_force, v=[.2, .4, .6])
        initial_momentum = np.array(self.system.analysis.linear_momentum())
        np.testing.assert_allclose(initial_momentum, np.copy(p.v) * p.mass,
                                   atol=mom_tol)
        while True:
            self.system.integrator.run(500)

            measured_momentum = self.system.analysis.linear_momentum()
            coupling_force = -(p.f - p.ext_force)
            compensation = -TIME_STEP / 2 * coupling_force

            np.testing.assert_allclose(measured_momentum + compensation,
                                       initial_momentum, atol=self.atol)
            if np.linalg.norm(p.f) < 0.01 \
               and np.all(np.abs(p.pos) > 10.1 * self.system.box_l):
                break

        # Make sure, the particle has crossed the periodic boundaries
        self.assertGreater(max(np.abs(p.v)) * self.system.time,
                           self.system.box_l[0])


@ut.skipIf(TestLBMomentumConservation.n_nodes == 1,
           "LB with regular decomposition already tested with 2 MPI ranks")
@utx.skipIfMissingFeatures(["WALBERLA", "EXTERNAL_FORCES"])
class TestLBMomentumConservationRegularWalberla(
        TestLBMomentumConservation, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": False}
    atol = 1.2e-4

    def set_cellsystem(self):
        self.system.cell_system.set_regular_decomposition()


@ut.skipIf(TestLBMomentumConservation.n_nodes == 1,
           "LB with regular decomposition already tested with 2 MPI ranks")
@utx.skipIfMissingFeatures(["WALBERLA", "EXTERNAL_FORCES"])
class TestLBMomentumConservationRegularWalberlaSinglePrecision(
        TestLBMomentumConservation, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": True}
    atol = 6.5e-4

    def set_cellsystem(self):
        self.system.cell_system.set_regular_decomposition()


@utx.skipIfMissingFeatures(["WALBERLA", "EXTERNAL_FORCES"])
class TestLBCPUMomentumConservationHybridNSquareWalberla(
        TestLBMomentumConservation, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": False}
    atol = 1.2e-4

    def set_cellsystem(self):
        self.system.cell_system.set_hybrid_decomposition(
            n_square_types={0}, cutoff_regular=1)


@utx.skipIfMissingFeatures(["WALBERLA", "EXTERNAL_FORCES"])
class TestLBCPUMomentumConservationHybridNSquareWalberlaSinglePrecision(
        TestLBMomentumConservation, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": True}
    atol = 6.5e-4

    def set_cellsystem(self):
        self.system.cell_system.set_hybrid_decomposition(
            n_square_types={0}, cutoff_regular=1)


@utx.skipIfMissingFeatures(["WALBERLA", "EXTERNAL_FORCES"])
class TestLBCPUMomentumConservationHybridRegularWalberla(
        TestLBMomentumConservation, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": False}
    atol = 1.2e-4

    def set_cellsystem(self):
        self.system.cell_system.set_hybrid_decomposition(
            n_square_types={1}, cutoff_regular=1)


@utx.skipIfMissingFeatures(["WALBERLA", "EXTERNAL_FORCES"])
class TestLBCPUMomentumConservationHybridRegularWalberlaSinglePrecision(
        TestLBMomentumConservation, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": True}
    atol = 6.5e-4

    def set_cellsystem(self):
        self.system.cell_system.set_hybrid_decomposition(
            n_square_types={1}, cutoff_regular=1)


@utx.skipIfMissingFeatures(["WALBERLA", "EXTERNAL_FORCES"])
class TestLBMomentumConservationNSquareWalberla(
        TestLBMomentumConservation, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": False}
    atol = 1.2e-4

    def set_cellsystem(self):
        self.system.cell_system.set_n_square()


@utx.skipIfMissingFeatures(["WALBERLA", "EXTERNAL_FORCES"])
class TestLBMomentumConservationNSquareWalberlaSinglePrecision(
        TestLBMomentumConservation, ut.TestCase):

    lb_class = espressomd.lb.LBFluidWalberla
    lb_params = {"single_precision": True}
    atol = 6.5e-4

    def set_cellsystem(self):
        self.system.cell_system.set_n_square()


if __name__ == "__main__":
    ut.main()
