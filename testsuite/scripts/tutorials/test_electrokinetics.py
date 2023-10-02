#
# Copyright (C) 2019-2023 The ESPResSo project
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
import importlib_wrapper as iw
import numpy as np
import scipy.optimize

tutorial, skipIfMissingFeatures = iw.configure_and_import(
    "@TUTORIALS_DIR@/electrokinetics/electrokinetics.py", TOTAL_FRAMES=0)


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    system = tutorial.system

    def test_analytic_solutions(self):
        # check advection-diffusion
        mu_simulated = tutorial.positions_diagonal[np.argmax(
            tutorial.values_diagonal)]
        mu_analytic = tutorial.positions_analytic[np.argmax(
            tutorial.values_analytic)]
        tol = 2. * tutorial.AGRID
        self.assertAlmostEqual(mu_simulated, mu_analytic, delta=tol)
        # check electroosmotic flow
        np.testing.assert_allclose(
            np.copy(tutorial.density_eof),
            tutorial.analytic_density_eof,
            rtol=0.005)
        np.testing.assert_allclose(
            np.copy(tutorial.velocity_eof),
            tutorial.analytic_velocity_eof,
            rtol=0.05)
        np.testing.assert_allclose(
            np.copy(tutorial.pressure_tensor_eof),
            tutorial.analytic_pressure_tensor_eof,
            rtol=0.05)
        # check pressure-driven flow
        np.testing.assert_allclose(
            np.copy(tutorial.density_pressure),
            tutorial.analytic_density_eof,
            rtol=0.005)
        np.testing.assert_allclose(
            np.copy(tutorial.velocity_pressure),
            tutorial.analytic_velocity_pressure,
            rtol=0.05)
        np.testing.assert_allclose(
            np.copy(tutorial.pressure_tensor_pressure),
            tutorial.analytic_pressure_tensor_pressure,
            rtol=0.05)

    def test_karman_vortex_street(self):
        """
        Check for the formation of a Kármán vortex street. Due to turbulence,
        a wavy pattern emerges.
        """
        def get_frequency_species(species):
            """
            Compute the principal wavevector of a chemical species.
            """
            fdata = np.zeros(16)
            for i in range(32, 80):
                rdata = species[i, :, 0].density
                fdata += np.abs(np.fft.fft(rdata - np.mean(rdata)))[:16]
            return np.argmax(fdata)

        def get_phase_karman(species):
            """
            Compute the time-dependent phase of a turbulent flow profile.
            """
            phase = []
            k = 2  # wavevector for product species
            for i in range(36, 68):
                rdata = species[i, :, 0].density
                phase.append(np.angle(np.fft.fft(rdata - np.mean(rdata))[k]))
            return np.array(phase)

        def cosine_kernel(x, magnitude, freq, phase):
            return magnitude * np.cos(x * freq + phase)

        tutorial.system.integrator.run(2000)
        # check for finite values
        for species in (*tutorial.educt_species, *tutorial.product_species):
            assert np.all(np.isfinite(species[:, :, :].density))
            assert np.all(species[:, :, :].density >= 0)
        assert np.all(np.isfinite(tutorial.lbf[:, :, :].velocity))
        # there is only one inlet per educt, thus wavelength == box width
        self.assertEqual(get_frequency_species(tutorial.educt_species[0]), 1)
        self.assertEqual(get_frequency_species(tutorial.educt_species[1]), 1)
        # reaction happens in the mixing region, thus the frequency doubles
        self.assertEqual(get_frequency_species(tutorial.product_species[0]), 2)
        # check for turbulence onset
        ref_params = np.array([2., 0.12, 1. / 2. * np.pi])
        sim_phase = get_phase_karman(tutorial.product_species[0])
        xdata = np.arange(sim_phase.shape[0])
        popt, _ = scipy.optimize.curve_fit(
            cosine_kernel, xdata, sim_phase, p0=ref_params,
            bounds=([-4., 0.08, 0.], [4., 0.24, 2. * np.pi]))
        fit_phase = cosine_kernel(xdata, *popt)
        rmsd = np.sqrt(np.mean(np.square(sim_phase - fit_phase)))
        self.assertAlmostEqual(popt[2], ref_params[2], delta=0.20)
        self.assertAlmostEqual(popt[1], ref_params[1], delta=0.02)
        self.assertAlmostEqual(popt[0], ref_params[0], delta=0.80)
        self.assertLess(rmsd / abs(popt[0]), 0.2)


if __name__ == "__main__":
    ut.main()
