#
# Copyright (C) 2019-2026 The ESPResSo project
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
    "@TUTORIALS_DIR@/electrokinetics/electrokinetics_part3_reactive_flow.py",
    TOTAL_FRAMES=0)


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    system = tutorial.system

    @classmethod
    def setUpClass(cls):
        # the animation was skipped (TOTAL_FRAMES=0), thus the flow field and
        # the species densities still have to be evolved
        tutorial.system.integrator.run(2000)

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
        ref_params = np.array([1., 0.08, 1. / 20. * np.pi])
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

    def test_product_profile(self):
        """
        No product can form as long as the two educt streams are separated.
        Further downstream, the vortex street mixes them and the product
        density grows.
        """
        profile = tutorial.mean_density_profile(tutorial.product_species[0])
        self.assertEqual(len(profile), tutorial.BOX_L[0])
        # essentially no product upstream of the obstacles
        self.assertLess(np.max(profile[:tutorial.BOX_L[0] // 10]),
                        1e-3 * np.max(profile))
        # the product density grows downstream
        self.assertLess(np.mean(profile[20:40]), np.mean(profile[60:80]))

    def test_reynolds_number(self):
        """
        The flow must be fast enough for the wake to become unsteady.
        """
        velocity_mean = np.nanmean(np.where(
            tutorial.boundary_mask, np.nan,
            tutorial.lbf[:, :, 0].velocity[..., 0]))
        reynolds_number = velocity_mean * 2. * tutorial.OBSTACLE_RADIUS / \
            tutorial.VISCOSITY_KINEMATIC
        self.assertGreater(reynolds_number, 5.)
        self.assertLess(reynolds_number, 50.)


if __name__ == "__main__":
    ut.main()
