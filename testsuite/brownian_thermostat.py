#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
# Tests particle property setters/getters
from __future__ import print_function
import unittest as ut
import espressomd
import numpy as np
from espressomd.interactions import FeneBond
from time import time
from espressomd.accumulators import Correlator
from espressomd.observables import ParticleVelocities, ParticleBodyAngularVelocities

#TODO: this test should be as close to langevin_thermostat.py as possible
#      these tests probably even could be merged
#      after an implementation of the fine inertial BD  
@ut.skipIf(espressomd.has_features("THERMOSTAT_IGNORE_NON_VIRTUAL") or
           not espressomd.has_features("BROWNIAN_DYNAMICS"),
           "Skipped because of the features set")
class BrownianThermostat(ut.TestCase):
    """Tests the velocity distribution created by the Brownian thermostat against
       the single component Maxwell distribution."""

    s = espressomd.System(box_l=[1.0, 1.0, 1.0])
    s.cell_system.set_n_square()
    s.cell_system.skin = 0.3
    s.seed = range(s.cell_system.get_state()["n_nodes"])
    if espressomd.has_features("PARTIAL_PERIODIC"):
        s.periodicity = 0,0,0


    @classmethod
    def setUpClass(cls):
        np.random.seed(42)

    def single_component_maxwell(self, x1, x2, kT):
        """Integrate the probability density from x1 to x2 using the trapez rule"""
        x = np.linspace(x1, x2, 1000)
        return np.trapz(np.exp(-x**2 / (2. * kT)), x) / \
            np.sqrt(2. * np.pi * kT)

    def check_velocity_distribution(self, vel, minmax, n_bins, error_tol, kT):
        """check the recorded particle distributions in vel againsta histogram with n_bins bins. Drop velocities outside minmax. Check individual histogram bins up to an accuracy of error_tol agaisnt the analytical result for kT."""
        for i in range(3):
            hist = np.histogram(
                vel[:, i], range=(-minmax, minmax), bins=n_bins, normed=False)
            data = hist[0] / float(vel.shape[0])
            bins = hist[1]
            for j in range(n_bins):
                found = data[j]
                expected = self.single_component_maxwell(
                    bins[j], bins[j + 1], kT)
                self.assertLessEqual(abs(found - expected), error_tol)

    def test_aa_verify_single_component_maxwell(self):
        """Verifies the normalization of the analytical expression."""
        self.assertLessEqual(
            abs(self.single_component_maxwell(-10, 10, 4.) - 1.), 1E-4)

    def test_global_brownian(self):
        """Test for global Brownian parameters."""
        N = 200
        s = self.s
        s.part.clear()
        s.time_step = 0.02
        
        # Place particles
        s.part.add(pos=np.random.random((N, 3)))
       
        # Enable rotation if compiled in
        if espressomd.has_features("ROTATION"): 
            s.part[:].rotation = 1,1,1

        kT = 2.3
        gamma = 1.5
        s.thermostat.set_brownian(kT=kT, gamma=gamma)
        
        # Warmup
        s.integrator.run(100)

        # Sampling
        loops = 30
        v_stored = np.zeros((N * loops, 3))
        omega_stored = np.zeros((N * loops, 3))
        for i in range(loops):
            s.integrator.run(2)
            v_stored[i * N:(i + 1) * N, :] = s.part[:].v
            if espressomd.has_features("ROTATION"):
                omega_stored[i * N:(i + 1) * N, :] = s.part[:].omega_body

        v_minmax = 5
        bins = 5
        error_tol = 0.015
        self.check_velocity_distribution(
            v_stored, v_minmax, bins, error_tol, kT)
        if espressomd.has_features("ROTATION"): 
            self.check_velocity_distribution(
                omega_stored, v_minmax, bins, error_tol, kT)

if __name__ == "__main__":
    ut.main()
