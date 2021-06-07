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
import numpy as np
import espressomd
import espressomd.reaction_ensemble


class ReactionEnsembleTest(ut.TestCase):

    """Test the core implementation of the constant pH reaction ensemble."""

    N0 = 40
    c0 = 0.00028
    type_HA = 0
    type_A = 1
    type_H = 5
    temperature = 1.0
    # choose target alpha not too far from 0.5 to get good statistics in a
    # small number of steps
    pKa_minus_pH = -0.2
    pH = 2
    pKa = pKa_minus_pH + pH
    Ka = 10**(-pKa)
    box_l = (N0 / c0)**(1.0 / 3.0)
    system = espressomd.System(box_l=[box_l, box_l, box_l])
    np.random.seed(69)  # make reaction code fully deterministic
    system.cell_system.skin = 0.4
    system.time_step = 0.01
    RE = espressomd.reaction_ensemble.ConstantpHEnsemble(
        temperature=1.0, exclusion_radius=1, seed=44, conversion_factor_from_mol_per_l_to_1_div_sigma_cubed=1.0 / 37.1, pH=pH)

    @classmethod
    def setUpClass(cls):
        cls.system.part.add(
            pos=np.random.random((2 * cls.N0, 3)) * cls.system.box_l,
            type=cls.N0 * [cls.type_A, cls.type_H])

        cls.RE.constant_pH = cls.pH
        cls.RE.set_volume(cls.box_l**3)

    @classmethod
    def tearDownClass(cls):
        cls.RE.delete_reaction(0)

    @classmethod
    def ideal_alpha(cls, pH):
        return 1.0 / (1 + 10**(cls.pKa - pH))

    def test_ideal_titration_curve(self):
        num_concurrent_reactions = 1
        RE = ReactionEnsembleTest.RE
        RE.add_reaction(
            gamma=ReactionEnsembleTest.Ka**num_concurrent_reactions,
            reactant_types=[ReactionEnsembleTest.type_HA],
            reactant_coefficients=[num_concurrent_reactions],
            product_types=[
                ReactionEnsembleTest.type_A,
                ReactionEnsembleTest.type_H],
            product_coefficients=[
                num_concurrent_reactions,
                num_concurrent_reactions],
            default_charges={
                ReactionEnsembleTest.type_HA: 0,
                ReactionEnsembleTest.type_A: -1,
                ReactionEnsembleTest.type_H: +1},
            product_index_protons=1)

        N0 = ReactionEnsembleTest.N0
        type_A = ReactionEnsembleTest.type_A
        type_H = ReactionEnsembleTest.type_H
        type_HA = ReactionEnsembleTest.type_HA
        system = ReactionEnsembleTest.system
        # chemical warmup - get close to chemical equilibrium before we start
        # sampling
        RE.reaction(40 * N0)

        average_NH = 0.0
        average_NHA = 0.0
        average_NA = 0.0
        num_samples = 1000
        for _ in range(num_samples):
            RE.reaction(10)
            average_NH += system.number_of_particles(type=type_H)
            average_NHA += system.number_of_particles(type=type_HA)
            average_NA += system.number_of_particles(type=type_A)
        average_NH /= float(num_samples)
        average_NA /= float(num_samples)
        average_NHA /= float(num_samples)
        average_alpha = average_NA / float(N0)
        # note you cannot calculate the pH via -log10(<NH>/volume) in the
        # constant pH ensemble, since the volume is totally arbitrary and does
        # not influence the average number of protons
        pH = ReactionEnsembleTest.pH
        pKa = ReactionEnsembleTest.pKa
        target_alpha = ReactionEnsembleTest.ideal_alpha(pH)
        rel_error_alpha = abs(average_alpha - target_alpha) / target_alpha
        self.assertLess(
            rel_error_alpha,
            0.015,
            msg="\nDeviation from ideal titration curve is too big for the given input parameters.\n"
            + f"  pH: {pH:.2f}"
            + f"  pKa: {pKa:.2f}"
            + f"  average NH: {average_NH:.1f}"
            + f"  average NA: {average_NA:.1f}"
            + f"  average NHA: {average_NHA:.1f}"
            + f"  average alpha: {average_alpha:.3f}"
            + f"  target alpha: {target_alpha:.3f}"
            + f"  rel_error: {rel_error_alpha * 100:.1f}%"
        )


if __name__ == "__main__":
    ut.main()
