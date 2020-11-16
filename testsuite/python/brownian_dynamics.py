#
# Copyright (C) 2013-2019 The ESPResSo project
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
from espressomd.accumulators import Correlator
from espressomd.observables import ParticlePositions
from thermostats_common import ThermostatsCommon


class ThermalizedBond(ut.TestCase, ThermostatsCommon):

    """Tests velocity distributions and diffusion for Brownian Dynamics"""
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    system.cell_system.set_domain_decomposition(use_verlet_lists=True)
    system.cell_system.skin = 0
    system.periodicity = [0, 0, 0]
    system.integrator.set_brownian_dynamics()

    def setUp(self):
        np.random.seed(42)

    def tearDown(self):
        self.system.time_step = 1e-12
        self.system.cell_system.skin = 0.0
        self.system.part.clear()
        self.system.auto_update_accumulators.clear()

    def check_vel_dist_global_temp(self, recalc_forces, loops):
        """Test velocity distribution for global temperature parameters.

        Parameters
        ----------
        recalc_forces : :obj:`bool`
            True if the forces should be recalculated after every step.
        loops : :obj:`int`
            Number of sampling loops
        """
        N = 200
        system = self.system
        system.time_step = 1.6
        kT = 1.1
        gamma = 3.5
        system.thermostat.set_brownian(kT=kT, gamma=gamma, seed=41)
        v_minmax = 5
        bins = 4
        error_tol = 0.01
        self.check_global(
            N, kT, loops, v_minmax, bins, error_tol, recalc_forces)

    def test_vel_dist_global_temp(self):
        """Test velocity distribution for global Brownian parameters."""
        self.check_vel_dist_global_temp(False, loops=200)

    def test_vel_dist_global_temp_initial_forces(self):
        """Test velocity distribution for global Brownian parameters,
           when using the initial force calculation.
        """
        self.check_vel_dist_global_temp(True, loops=170)

    @utx.skipIfMissingFeatures("BROWNIAN_PER_PARTICLE")
    def test_vel_dist_per_particle(self):
        """Test Brownian dynamics with particle-specific kT and gamma. Covers
           all combinations of particle-specific gamma and temp set or not set.
        """
        N = 400
        system = self.system
        system.time_step = 1.9
        kT = 0.9
        gamma = 3.2
        gamma2 = 4.3
        kT2 = 1.5
        system.thermostat.set_brownian(kT=kT, gamma=gamma, seed=41)
        loops = 200
        v_minmax = 5
        bins = 4
        error_tol = 0.012
        self.check_per_particle(
            N, kT, kT2, gamma2, loops, v_minmax, bins, error_tol)

    def setup_diff_mass_rinertia(self, p):
        if espressomd.has_features("MASS"):
            p.mass = 0.5
        if espressomd.has_features("ROTATION"):
            p.rotation = [1, 1, 1]
            # Make sure rinertia does not change diff coeff
            if espressomd.has_features("ROTATIONAL_INERTIA"):
                p.rinertia = [0.4, 0.4, 0.4]

    def test_msd_global_temp(self):
        """Tests diffusion via MSD for global gamma and temperature"""

        gamma = 9.4
        kT = 0.37
        dt = 0.5

        system = self.system
        p = system.part.add(pos=(0, 0, 0), id=0)
        system.time_step = dt
        system.thermostat.set_brownian(kT=kT, gamma=gamma, seed=41)
        system.cell_system.skin = 0.4

        pos_obs = ParticlePositions(ids=(p.id,))

        c_pos = Correlator(obs1=pos_obs, tau_lin=16, tau_max=100., delta_N=1,
                           corr_operation="square_distance_componentwise",
                           compress1="discard1")
        system.auto_update_accumulators.add(c_pos)

        system.integrator.run(30000)

        c_pos.finalize()

        # Check MSD
        msd = c_pos.result()
        tau = c_pos.lag_times()
        system.auto_update_accumulators.clear()

        def expected_msd(x):
            return 2. * kT / gamma * x

        for i in range(2, 6):
            np.testing.assert_allclose(
                msd[i], expected_msd(tau[i]), rtol=0.02)

    @utx.skipIfMissingFeatures("VIRTUAL_SITES")
    def test_07__virtual(self):
        system = self.system
        system.time_step = 0.01

        virtual = system.part.add(pos=[0, 0, 0], virtual=True, v=[1, 0, 0])
        physical = system.part.add(pos=[0, 0, 0], virtual=False, v=[1, 0, 0])

        system.thermostat.set_brownian(
            kT=0, gamma=1, gamma_rotation=1., act_on_virtual=False, seed=41)

        system.integrator.run(1)

        np.testing.assert_almost_equal(np.copy(virtual.v), [1, 0, 0])
        np.testing.assert_almost_equal(np.copy(physical.v), [0, 0, 0])

        virtual.pos = physical.pos = [0, 0, 0]
        system.thermostat.set_brownian(
            kT=0, gamma=1, gamma_rotation=1., act_on_virtual=True, seed=41)
        system.integrator.run(1)

        np.testing.assert_almost_equal(np.copy(virtual.v), [0, 0, 0])
        np.testing.assert_almost_equal(np.copy(physical.v), [0, 0, 0])

    def test_08__noise_correlation(self):
        """Checks that the Brownian noise is uncorrelated"""

        system = self.system
        system.time_step = 0.01
        system.cell_system.skin = 0.1
        kT = 3.2
        system.thermostat.set_brownian(kT=kT, gamma=2.1, seed=17)
        system.part.add(id=(0, 1), pos=np.zeros((2, 3)))
        steps = int(1e4)
        error_delta = 0.04
        self.check_noise_correlation(kT, steps, error_delta)


if __name__ == "__main__":
    ut.main()
