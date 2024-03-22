#
# Copyright (C) 2023 The ESPResSo project
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

import espressomd
import espressomd.reaction_methods
import unittest as ut
import numpy as np


class Test(ut.TestCase):

    """Test the trivial reaction 1A <-> 1B"""

    system = espressomd.System(box_l=[1., 1., 1.])
    system.time_step = 0.1
    system.cell_system.skin = 0.1
    np.random.seed(42)

    def test_equilibration(self):
        system = self.system
        N0 = 40
        c0 = 0.001
        types = {"A": 0, "B": 1}
        system.box_l = np.ones(3) * np.cbrt(N0 / c0)
        RE = espressomd.reaction_methods.ReactionEnsemble(
            seed=42, kT=1., exclusion_range=1., search_algorithm="parallel")
        RE.set_non_interacting_type(type=max(types.values()) + 1)
        system.part.add(
            pos=np.random.random((N0, 3)) * system.box_l,
            type=N0 * [types["A"]])
        RE.add_reaction(
            gamma=1.,
            reactant_types=[types["A"]],
            reactant_coefficients=[1],
            product_types=[types["B"]],
            product_coefficients=[1],
            default_charges={types["A"]: 0, types["B"]: 0})

        # warmup
        system.thermostat.set_langevin(kT=1., gamma=3., seed=42)
        system.integrator.run(200)
        RE.reaction(steps=5 * N0)

        # sampling
        average_NA = 0.0
        num_samples = 100
        for _ in range(num_samples):
            RE.reaction(steps=10)
            system.integrator.run(20)
            average_NA += system.number_of_particles(type=types["A"])
        average_NA /= num_samples

        alpha = average_NA / float(N0)
        rate0 = RE.get_acceptance_rate_reaction(reaction_id=0)
        rate1 = RE.get_acceptance_rate_reaction(reaction_id=1)
        self.assertAlmostEqual(alpha, 0.50, delta=0.01)
        self.assertAlmostEqual(rate0, 0.85, delta=0.05)
        self.assertAlmostEqual(rate1, 0.85, delta=0.05)


if __name__ == "__main__":
    ut.main()
