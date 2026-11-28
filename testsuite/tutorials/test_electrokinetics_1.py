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

tutorial, skipIfMissingFeatures = iw.configure_and_import(
    "@TUTORIALS_DIR@/electrokinetics/electrokinetics_part1_advection_diffusion.py")


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    system = tutorial.system

    def test_droplet_position(self):
        """
        The peak of the advected droplet must sit where the fundamental
        solution of the advection-diffusion equation predicts it.
        """
        mu_simulated = tutorial.positions_diagonal[np.argmax(
            tutorial.values_diagonal)]
        mu_analytic = tutorial.positions_analytic[np.argmax(
            tutorial.values_analytic)]
        tol = 2. * tutorial.AGRID
        self.assertAlmostEqual(mu_simulated, mu_analytic, delta=tol)

    def test_conservation_and_drift(self):
        """
        The finite-difference scheme is flux-based and therefore conserves
        the total amount of solute exactly; the center of mass is advected
        with the flow velocity.
        """
        self.assertAlmostEqual(tutorial.total_amount, 1., delta=1e-8)
        np.testing.assert_allclose(np.copy(tutorial.com), tutorial.drift,
                                   atol=tutorial.AGRID)

    def test_numerical_diffusion(self):
        """
        The discretization of the advection term adds an artificial
        contribution to the diffusion coefficient.
        """
        self.assertGreater(tutorial.diffusion_effective,
                           tutorial.DIFFUSION_COEFFICIENT)
        self.assertLess(tutorial.diffusion_effective,
                        2. * tutorial.DIFFUSION_COEFFICIENT)


if __name__ == "__main__":
    ut.main()
