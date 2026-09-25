#
# Copyright (C) 2013-2026 The ESPResSo project
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
import espressomd.reaction_methods
import numpy as np


@utx.skipIfMissingFeatures(["ELECTROSTATICS"])
class ReactionMethodsBookkeepingTest(ut.TestCase):
    """
    Test that two different instances of the reaction methods
    do not break the particle id bookkeeping.
    """

    def test_reaction_bookeeping(self):
        pH = 10.
        pKa = 7.
        exclusion_range = 1.
        seed = 12345
        kT = 1.
        BOX_LENGTH = 100.
        N_SALT = 10
        N_acid = 10

        types = {"H": 0, "Na": 1, "Cl": 2, "HA": 3, "A": 4}
        charges = {"H": 1.0, "Na": 1.0, "Cl": -1.0, "HA": 0.0, "A": -1.0}
        system = espressomd.System(box_l=[BOX_LENGTH, ] * 3)

        cph = espressomd.reaction_methods.ConstantpHEnsemble(
            system=system,
            constant_pH=pH,
            kT=kT,
            exclusion_range=exclusion_range,
            seed=seed)
        widom = espressomd.reaction_methods.WidomInsertion(
            system=system,
            kT=kT,
            seed=seed)
        system.part.add(type=[types["Na"]] * N_SALT,
                        pos=np.random.rand(N_SALT, 3) * BOX_LENGTH,
                        q=[charges["Na"]] * N_SALT,
                        id=list(range(20, 20 + N_SALT)))
        system.part.add(type=[types["Cl"]] * N_SALT,
                        pos=np.random.rand(N_SALT, 3) * BOX_LENGTH,
                        q=[charges["Cl"]] * N_SALT)
        system.part.add(type=[types["HA"]] * N_acid,
                        pos=np.random.rand(N_acid, 3) * BOX_LENGTH,
                        q=[charges["HA"]] * N_acid)

        cph.add_reaction(
            gamma=10**(-pKa),
            reactant_types=[types["HA"]],
            product_types=[types["A"], types["H"]],
            default_charges={types["HA"]: charges["HA"],
                             types["A"]: charges["A"],
                             types["H"]: charges["H"]}
        )

        widom.add_reaction(
            reactant_types=[],
            reactant_coefficients=[],
            product_types=[types["Na"], types["Cl"]],
            product_coefficients=[1, 1],
            default_charges={types["Na"]: charges["Na"],
                             types["Cl"]: charges["Cl"]}
        )

        widom.calculate_particle_insertion_potential_energy(reaction_id=0)
        cph.reaction(steps=100)

        # Measure the chemical potential
        for _ in range(50):
            widom.calculate_particle_insertion_potential_energy(reaction_id=0)
            total_charge = sum(system.part.all().q)
            self.assertEqual(total_charge, 0.)


if __name__ == "__main__":
    ut.main()
