#
# Copyright (C) 2013-2021 The ESPResSo project
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
import numpy as np
import espressomd
import espressomd.reaction_ensemble


class ConstantpHTest(ut.TestCase):

    """Test the core implementation of the constant pH reaction ensemble."""

    np.random.seed(42)

    def test_ideal_alpha(self):
        N0 = 40
        c0 = 0.00028
        type_HA = 0
        type_A = 1
        type_H = 5
        pH = 2.0
        pKa = 2.5
        box_l = (N0 / c0)**(1.0 / 3.0)
        system = espressomd.System(box_l=3 * [box_l])
        system.cell_system.skin = 0.4
        system.time_step = 0.01
        system.part.add(pos=np.random.random((2 * N0, 3)) * system.box_l,
                        type=N0 * [type_A, type_H])

        RE = espressomd.reaction_ensemble.ConstantpHEnsemble(
            temperature=1.0,
            exclusion_radius=1,
            seed=44)
        RE.add_reaction(
            gamma=10**(-pKa),
            reactant_types=[type_HA],
            reactant_coefficients=[1],
            product_types=[type_A, type_H],
            product_coefficients=[1, 1],
            default_charges={type_HA: 0, type_A: -1, type_H: +1})
        RE.constant_pH = pH

        # equilibration
        RE.reaction(800)

        # sampling
        alphas = []
        for _ in range(80):
            RE.reaction(15)
            num_H = system.number_of_particles(type=type_H)
            num_HA = system.number_of_particles(type=type_HA)
            num_A = system.number_of_particles(type=type_A)
            self.assertEqual(num_H, num_A)
            self.assertEqual(num_A + num_HA, N0)
            alphas.append(num_A / N0)

        # note: you cannot calculate the pH via -log10(<NH>/volume) in the
        # constant pH ensemble, since the volume is totally arbitrary and does
        # not influence the average number of protons
        target_alpha = 1. / (1. + 10.**(pKa - pH))
        alpha_mean = np.mean(alphas)
        alpha_error = np.std(alphas) / np.sqrt(len(alphas))
        self.assertAlmostEqual(alpha_mean, target_alpha, delta=0.015)
        self.assertAlmostEqual(alpha_error, 0., delta=0.01)


if __name__ == "__main__":
    ut.main()
