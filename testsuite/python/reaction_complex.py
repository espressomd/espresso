#
# Copyright (C) 2013-2022 The ESPResSo project
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
import scipy.optimize
import espressomd
import espressomd.reaction_methods


class Test(ut.TestCase):

    """
    Test the reaction ensemble for reaction 2A + 3B <-> 4C + 1D + 3E.
    """

    system = espressomd.System(box_l=3 * [35.])
    system.time_step = 0.02
    system.cell_system.skin = 0.4
    np.random.seed(seed=42)

    def test_equilibrium_concentrations(self):
        # reactant and product types
        type_A = 0
        type_B = 1
        type_C = 2
        type_D = 3
        type_E = 4
        types = [type_A, type_B, type_C, type_D, type_E]
        # reactant and product stoichiometric coefficients
        nu_A = -2
        nu_B = -3
        nu_C = 4
        nu_D = 1
        nu_E = 3
        nu_s = [nu_A, nu_B, nu_C, nu_D, nu_E]
        # reaction constant
        K = 0.001
        # reference concentration: 1 mol/l
        c_ref = 1.0
        # simulation units: 1 sigma = 3.55 Angstrom
        conv_sim_unit_to_mol_per_l = 37.1
        Gamma = K * (c_ref / conv_sim_unit_to_mol_per_l)**np.sum(nu_s)
        volume = self.system.volume()

        # setup N0 batches of reactants
        N0 = 30
        self.system.part.add(
            pos=np.random.random((-nu_A * N0, 3)) * self.system.box_l,
            type=-nu_A * N0 * [type_A])
        self.system.part.add(
            pos=np.random.random((-nu_B * N0, 3)) * self.system.box_l,
            type=-nu_B * N0 * [type_B])

        # use an exclusion radius of 0 to simulate an ideal gas
        RE = espressomd.reaction_methods.ReactionEnsemble(
            kT=1., exclusion_range=0., seed=42)

        RE.add_reaction(
            gamma=Gamma, reactant_types=[type_A, type_B],
            reactant_coefficients=[abs(nu_A), abs(nu_B)],
            product_types=[type_C, type_D, type_E],
            product_coefficients=[nu_C, nu_D, nu_E],
            default_charges={type_A: 0, type_B: 0, type_C: 0, type_D: 0, type_E: 0})

        # Set the hidden particle type to the lowest possible number to speed
        # up the simulation
        RE.set_non_interacting_type(type=max(types) + 1)

        # warmup
        RE.reaction(reaction_steps=50)

        numbers = {type_A: [], type_B: [], type_C: [], type_D: [], type_E: []}
        for _ in range(40):
            RE.reaction(reaction_steps=6)
            for key in types:
                numbers[key].append(self.system.number_of_particles(type=key))

        def calculate_K(c_s, nu_s):
            return np.prod(
                [(c_i / c_ref)**nu_i for c_i, nu_i in zip(c_s, nu_s)])

        def equations(concentrations):
            c_A, c_B, c_C, c_D, c_E = concentrations
            eq1 = K - calculate_K(concentrations, nu_s)
            eq2 = N0 - (1.0 / abs(nu_A) * c_A / conv_sim_unit_to_mol_per_l +
                        1.0 / nu_D * c_D / conv_sim_unit_to_mol_per_l) * volume
            eq3 = c_A / c_B - float(nu_A) / float(nu_B)
            eq4 = c_C / c_D - float(nu_C) / float(nu_D)
            eq5 = c_C / c_E - float(nu_C) / float(nu_E)
            return (eq1, eq2, eq3, eq4, eq5)

        concentrations_sim = {
            type_i: np.mean(n_i) / volume * conv_sim_unit_to_mol_per_l
            for type_i, n_i in numbers.items()}
        concentrations_analytic = dict(zip(types, scipy.optimize.fsolve(
            equations, list(concentrations_sim.values()))))

        K_sim = calculate_K(concentrations_sim.values(), nu_s)
        N0_sim = (1.0 / abs(nu_A) * concentrations_sim[type_A] + 1.0 / nu_D *
                  concentrations_sim[type_D]) / conv_sim_unit_to_mol_per_l * volume

        err_msg = "Concentration of species {} doesn't match analytical result"
        for key in types:
            self.assertAlmostEqual(concentrations_sim[key],
                                   concentrations_analytic[key],
                                   delta=1e-3,
                                   msg=err_msg.format(key))
        self.assertAlmostEqual(K_sim, K, delta=1e-3)
        self.assertAlmostEqual(N0_sim, N0, delta=1e-3)


if __name__ == "__main__":
    ut.main()
