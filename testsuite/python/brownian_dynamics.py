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
from espressomd.accumulators import Correlator, TimeSeries
from espressomd.observables import ParticleVelocities, ParticleBodyAngularVelocities, ParticlePositions
from tests_common import single_component_maxwell


class BrownianDynamics(ut.TestCase):

    """Tests velocity distributions and diffusion for Brownian Dynamics"""
    system = espressomd.System(box_l=[1.0, 1.0, 1.0])
    system.cell_system.set_domain_decomposition(use_verlet_lists=True)
    system.cell_system.skin = 0
    system.periodicity = [0, 0, 0]
    system.integrator.set_brownian_dynamics()

    @classmethod
    def setUpClass(cls):
        np.random.seed(42)

    def check_velocity_distribution(self, vel, minmax, n_bins, error_tol, kT):
        """check the recorded particle distributions in velocity against a
           histogram with n_bins bins. Drop velocities outside minmax. Check
           individual histogram bins up to an accuracy of error_tol against the
           analytical result for kT."""
        for i in range(3):
            hist = np.histogram(
                vel[:, i], range=(-minmax, minmax), bins=n_bins, density=False)
            data = hist[0] / float(vel.shape[0])
            bins = hist[1]
            for j in range(n_bins):
                found = data[j]
                expected = single_component_maxwell(bins[j], bins[j + 1], kT)
                self.assertLessEqual(abs(found - expected), error_tol)

    def test_00_verify_single_component_maxwell(self):
        """Verifies the normalization of the analytical expression."""
        self.assertLessEqual(
            abs(single_component_maxwell(-10, 10, 4.) - 1.), 1E-4)

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
        system.part.clear()
        system.time_step = 1.6

        # Place particles
        system.part.add(pos=np.random.random((N, 3)))

        # Enable rotation if compiled in
        if espressomd.has_features("ROTATION"):
            system.part[:].rotation = [1, 1, 1]

        kT = 1.1
        gamma = 3.5
        system.thermostat.set_brownian(kT=kT, gamma=gamma, seed=41)

        # Warmup
        system.integrator.run(20)

        # Sampling
        v_stored = np.zeros((N * loops, 3))
        omega_stored = np.zeros((N * loops, 3))
        for i in range(loops):
            system.integrator.run(1, recalc_forces=recalc_forces)
            v_stored[i * N:(i + 1) * N, :] = system.part[:].v
            if espressomd.has_features("ROTATION"):
                omega_stored[i * N:(i + 1) * N, :] = system.part[:].omega_body

        v_minmax = 5
        bins = 4
        error_tol = 0.01
        self.check_velocity_distribution(
            v_stored, v_minmax, bins, error_tol, kT)
        if espressomd.has_features("ROTATION"):
            self.check_velocity_distribution(
                omega_stored, v_minmax, bins, error_tol, kT)

    def test_vel_dist_global_temp(self):
        """Test velocity distribution for global temperature."""
        self.check_vel_dist_global_temp(False, loops=1500)

    def test_vel_dist_global_temp_initial_forces(self):
        """Test velocity distribution for global Brownian parameters,
           when using the initial force calculation.
        """
        self.check_vel_dist_global_temp(True, loops=170)

    @utx.skipIfMissingFeatures("BROWNIAN_PER_PARTICLE")
    def test_05__brownian_per_particle(self):
        """Test Brownian dynamics with particle specific kT and gamma. Covers all combinations of
           particle specific gamma and temp set or not set.
        """
        N = 400
        system = self.system
        system.part.clear()
        system.time_step = 1.9
        system.part.add(pos=np.random.random((N, 3)))
        if espressomd.has_features("ROTATION"):
            system.part[:].rotation = [1, 1, 1]

        kT = 0.9
        gamma = 3.2
        gamma2 = 14.3
        kT2 = 1.5
        system.thermostat.set_brownian(kT=kT, gamma=gamma, seed=41)
        # Set different kT on 2nd half of particles
        system.part[int(N / 2):].temp = kT2
        # Set different gamma on half of the particles (overlap over both kTs)
        if espressomd.has_features("PARTICLE_ANISOTROPY"):
            system.part[int(N / 4):int(3 * N / 4)].gamma = 3 * [gamma2]
        else:
            system.part[int(N / 4):int(3 * N / 4)].gamma = gamma2

        system.integrator.run(50)
        loops = 300

        v_kT = np.zeros((int(N / 2) * loops, 3))
        v_kT2 = np.zeros((int(N / 2 * loops), 3))

        if espressomd.has_features("ROTATION"):
            omega_kT = np.zeros((int(N / 2) * loops, 3))
            omega_kT2 = np.zeros((int(N / 2 * loops), 3))

        for i in range(loops):
            system.integrator.run(1)
            v_kT[int(i * N / 2):int((i + 1) * N / 2),
                 :] = system.part[:int(N / 2)].v
            v_kT2[int(i * N / 2):int((i + 1) * N / 2),
                  :] = system.part[int(N / 2):].v

            if espressomd.has_features("ROTATION"):
                omega_kT[int(i * N / 2):int((i + 1) * N / 2), :] = \
                    system.part[:int(N / 2)].omega_body
                omega_kT2[int(i * N / 2):int((i + 1) * N / 2), :] = \
                    system.part[int(N / 2):].omega_body
        v_minmax = 5
        bins = 4
        error_tol = 0.012
        self.check_velocity_distribution(v_kT, v_minmax, bins, error_tol, kT)
        self.check_velocity_distribution(v_kT2, v_minmax, bins, error_tol, kT2)

        if espressomd.has_features("ROTATION"):
            self.check_velocity_distribution(
                omega_kT, v_minmax, bins, error_tol, kT)
            self.check_velocity_distribution(
                omega_kT2, v_minmax, bins, error_tol, kT2)

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
        system.part.clear()
        p = system.part.add(pos=(0, 0, 0), id=0)
        system.time_step = dt
        system.thermostat.set_brownian(kT=kT, gamma=gamma, seed=42)
        system.cell_system.skin = 0.4

        pos_obs = ParticlePositions(ids=(p.id,))

        c_pos = Correlator(obs1=pos_obs, tau_lin=16, tau_max=100., delta_N=10,
                           corr_operation="square_distance_componentwise",
                           compress1="discard1")
        system.auto_update_accumulators.add(c_pos)

        system.integrator.run(500000)

        c_pos.finalize()

        # Check MSD
        msd = c_pos.result()
        system.auto_update_accumulators.clear()

        def expected_msd(x):
            return 2. * kT / gamma * x

        for i in range(2, 6):
            np.testing.assert_allclose(
                msd[i, 2:5], expected_msd(msd[i, 0]), rtol=0.02)

    @utx.skipIfMissingFeatures("VIRTUAL_SITES")
    def test_07__virtual(self):
        system = self.system
        system.time_step = 0.01
        system.part.clear()

        virtual = system.part.add(pos=[0, 0, 0], virtual=True, v=[1, 0, 0])
        physical = system.part.add(pos=[0, 0, 0], virtual=False, v=[1, 0, 0])

        system.thermostat.set_brownian(
            kT=0, gamma=1, gamma_rotation=1., act_on_virtual=False, seed=41)

        system.integrator.run(1)

        np.testing.assert_almost_equal(np.copy(virtual.v), [1, 0, 0])
        np.testing.assert_almost_equal(np.copy(physical.v), [0, 0, 0])

        system.part.clear()
        virtual = system.part.add(pos=[0, 0, 0], virtual=True, v=[1, 0, 0])
        physical = system.part.add(pos=[0, 0, 0], virtual=False, v=[1, 0, 0])

        system.thermostat.set_brownian(
            kT=0, gamma=1, gamma_rotation=1., act_on_virtual=True, seed=41)
        system.integrator.run(1)

        np.testing.assert_almost_equal(np.copy(virtual.v), [0, 0, 0])
        np.testing.assert_almost_equal(np.copy(physical.v), [0, 0, 0])

    def test_08__noise_correlation(self):
        """Checks that the Brownian noise is uncorrelated"""

        system = self.system
        system.part.clear()
        system.time_step = 0.01
        system.cell_system.skin = 0.1
        system.part.add(id=(0, 1), pos=np.zeros((2, 3)))
        vel_obs = ParticleVelocities(ids=system.part[:].id)
        vel_series = TimeSeries(obs=vel_obs)
        system.auto_update_accumulators.add(vel_series)
        if espressomd.has_features("ROTATION"):
            system.part[:].rotation = (1, 1, 1)
            omega_obs = ParticleBodyAngularVelocities(ids=system.part[:].id)
            omega_series = TimeSeries(obs=omega_obs)
            system.auto_update_accumulators.add(omega_series)

        kT = 3.2
        system.thermostat.set_brownian(kT=kT, gamma=2.1, seed=17)
        steps = int(1e6)
        system.integrator.run(steps)
        system.auto_update_accumulators.clear()

        # test translational noise correlation
        vel = np.array(vel_series.time_series())
        for ind in range(2):
            for i in range(3):
                for j in range(i, 3):
                    corrcoef = np.dot(
                        vel[:, ind, i], vel[:, ind, j]) / steps / kT
                    if i == j:
                        self.assertAlmostEqual(corrcoef, 1.0, delta=0.04)
                    else:
                        self.assertLessEqual(np.abs(corrcoef), 0.04)

        # test rotational noise correlation
        if espressomd.has_features("ROTATION"):
            omega = np.array(omega_series.time_series())
            for ind in range(2):
                for i in range(3):
                    for j in range(3):
                        corrcoef = np.dot(
                            omega[:, ind, i], omega[:, ind, j]) / steps / kT
                        if i == j:
                            self.assertAlmostEqual(corrcoef, 1.0, delta=0.04)
                        else:
                            self.assertLessEqual(np.abs(corrcoef), 0.04)
                        # translational and angular velocities should be
                        # independent
                        corrcoef = np.dot(
                            vel[:, ind, i], omega[:, ind, j]) / steps / kT
                        self.assertLessEqual(np.abs(corrcoef), 0.04)


if __name__ == "__main__":
    ut.main()
