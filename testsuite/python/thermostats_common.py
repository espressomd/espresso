#
# Copyright (C) 2013-2020 The ESPResSo project
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

import espressomd
import espressomd.accumulators
import espressomd.observables


def single_component_maxwell(x1, x2, kT):
    """Integrate the probability density from x1 to x2 using the trapezoidal rule"""
    x = np.linspace(x1, x2, 1000)
    return np.trapz(np.exp(-x**2 / (2. * kT)), x) / np.sqrt(2. * np.pi * kT)


class ThermostatsCommon:

    """Tests the velocity distribution created by a thermostat."""

    def check_velocity_distribution(self, vel, minmax, n_bins, error_tol, kT):
        """Check the recorded particle distributions in velocity against a
           histogram with n_bins bins. Drop velocities outside minmax. Check
           individual histogram bins up to an accuracy of error_tol against
           the analytical result for kT."""
        for i in range(3):
            hist = np.histogram(
                vel[:, i], range=(-minmax, minmax), bins=n_bins, density=False)
            data = hist[0] / float(vel.shape[0])
            bins = hist[1]
            expected = []
            for j in range(n_bins):
                expected.append(single_component_maxwell(
                    bins[j], bins[j + 1], kT))
            np.testing.assert_allclose(data[:n_bins], expected, atol=error_tol)

    def test_00_verify_single_component_maxwell(self):
        """Verifies the normalization of the analytical expression."""
        self.assertAlmostEqual(
            single_component_maxwell(-10, 10, 4.), 1., delta=1E-4)

    def check_global(self, N, kT, steps, v_minmax,
                     n_bins, error_tol, recalc_forces):
        """
        Test Langevin/Brownian dynamics velocity distribution with global kT
        and gamma.

        Parameters
        ----------
        N : :obj:`int`
            Number of particles.
        kT : :obj:`float`
            Global temperature.
        steps : :obj:`int`
            Number of sampling steps.
        v_minmax : :obj:`float`
            Velocity range.
        n_bins : :obj:`int`
            Number of bins.
        error_tol : :obj:`float`
            Error tolerance.
        recalc_forces : :obj:`bool`
            True if the forces should be recalculated after every step.
        """
        system = self.system

        # Place particles
        partcls = system.part.add(pos=np.random.random((N, 3)))

        # Enable rotation if compiled in
        if espressomd.has_features("ROTATION"):
            partcls.rotation = 3 * [True]

        # Warmup
        system.integrator.run(20)

        vel_obs = espressomd.observables.ParticleVelocities(
            ids=partcls.id)
        vel_acc = espressomd.accumulators.TimeSeries(obs=vel_obs)
        system.auto_update_accumulators.add(vel_acc)

        if espressomd.has_features("ROTATION"):
            omega_obs = espressomd.observables.ParticleBodyAngularVelocities(
                ids=partcls.id)
            omega_acc = espressomd.accumulators.TimeSeries(obs=omega_obs)
            system.auto_update_accumulators.add(omega_acc)

        # Sampling
        v_stored = np.zeros((steps, N, 3))
        omega_stored = np.zeros((steps, N, 3))
        for i in range(steps):
            system.integrator.run(1, recalc_forces=recalc_forces)
            v_stored[i] = partcls.v
            if espressomd.has_features("ROTATION"):
                omega_stored[i] = partcls.omega_body

        vel = vel_acc.time_series().reshape((-1, 3))
        self.check_velocity_distribution(vel, v_minmax, n_bins, error_tol, kT)
        if espressomd.has_features("ROTATION"):
            omega = omega_acc.time_series().reshape((-1, 3))
            self.check_velocity_distribution(
                omega, v_minmax, n_bins, error_tol, kT)

    def check_per_particle(self, N, kT, gamma_local, steps, v_minmax, n_bins,
                           error_tol):
        """
        Test Langevin/Brownian dynamics velocity distribution with global and
        particle-specific kT/gamma. Covers all combinations of particle
        specific-gamma and temp set or not set.

        Parameters
        ----------
        N : :obj:`int`
            Number of particles.
        kT : :obj:`float`
            Global temperature.
        gamma_local : :obj:`float`
            Per-particle gamma.
        steps : :obj:`int`
            Number of sampling steps.
        v_minmax : :obj:`float`
            Velocity range.
        n_bins : :obj:`int`
            Number of bins.
        error_tol : :obj:`float`
            Error tolerance.
        """
        system = self.system
        partcls = system.part.add(pos=np.random.random((N, 3)))
        if espressomd.has_features("ROTATION"):
            partcls.rotation = [1, 1, 1]

        if espressomd.has_features("PARTICLE_ANISOTROPY"):
            gamma_local = 3 * [gamma_local]

        # Set different gamma on 2nd half of particles
        system.part.by_ids(range(N // 2, N)).gamma = gamma_local

        system.integrator.run(50)

        vel_obs = espressomd.observables.ParticleVelocities(
            ids=partcls.id)
        vel_acc = espressomd.accumulators.TimeSeries(obs=vel_obs)
        system.auto_update_accumulators.add(vel_acc)

        if espressomd.has_features("ROTATION"):
            omega_obs = espressomd.observables.ParticleBodyAngularVelocities(
                ids=partcls.id)
            omega_acc = espressomd.accumulators.TimeSeries(obs=omega_obs)
            system.auto_update_accumulators.add(omega_acc)

        system.integrator.run(steps)

        vel = vel_acc.time_series().reshape((-1, 3))
        self.check_velocity_distribution(vel, v_minmax, n_bins, error_tol, kT)
        if espressomd.has_features("ROTATION"):
            omega = omega_acc.time_series().reshape((-1, 3))
            self.check_velocity_distribution(
                omega, v_minmax, n_bins, error_tol, kT)

    def check_noise_correlation(self, kT, steps, delta):
        """Test the Langevin/Brownian noise is uncorrelated.

        Parameters
        ----------
        kT : :obj:`float`
            Global temperature.
        steps : :obj:`int`
            Number of sampling steps.
        delta : :obj:`float`
            Error tolerance.
        """

        system = self.system
        partcls = system.part.all()
        vel_obs = espressomd.observables.ParticleVelocities(
            ids=partcls.id)
        vel_series = espressomd.accumulators.TimeSeries(obs=vel_obs)
        system.auto_update_accumulators.add(vel_series)
        if espressomd.has_features("ROTATION"):
            partcls.rotation = 3 * [True]
            omega_obs = espressomd.observables.ParticleBodyAngularVelocities(
                ids=partcls.id)
            omega_series = espressomd.accumulators.TimeSeries(obs=omega_obs)
            system.auto_update_accumulators.add(omega_series)

        system.integrator.run(steps)

        # test translational noise correlation
        vel = np.array(vel_series.time_series())
        for ind in range(2):
            for i in range(3):
                for j in range(i, 3):
                    corrcoef = np.dot(
                        vel[:, ind, i], vel[:, ind, j]) / steps / kT
                    if i == j:
                        self.assertAlmostEqual(corrcoef, 1.0, delta=delta)
                    else:
                        self.assertAlmostEqual(corrcoef, 0.0, delta=delta)

        # test rotational noise correlation
        if espressomd.has_features("ROTATION"):
            omega = np.array(omega_series.time_series())
            for ind in range(2):
                for i in range(3):
                    for j in range(3):
                        corrcoef = np.dot(
                            omega[:, ind, i], omega[:, ind, j]) / steps / kT
                        if i == j:
                            self.assertAlmostEqual(corrcoef, 1.0, delta=delta)
                        else:
                            self.assertAlmostEqual(corrcoef, 0.0, delta=delta)
                        # translational and angular velocities should be
                        # independent
                        corrcoef = np.dot(
                            vel[:, ind, i], omega[:, ind, j]) / steps / kT
                        self.assertAlmostEqual(corrcoef, 0.0, delta=delta)
