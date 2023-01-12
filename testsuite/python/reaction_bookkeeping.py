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
        constant_pH=pH,
        kT=kT,
        exclusion_range=exclusion_range,
        seed=seed)
    widom = espressomd.reaction_methods.WidomInsertion(
        kT=kT,
        seed=seed)

    @classmethod
    def setUpClass(cls):
        cls.system.part.add(type=[cls.types["Na"]] * cls.N_SALT,
                            pos=np.random.rand(cls.N_SALT, 3) * cls.BOX_LENGTH,
                            q=[cls.charges["Na"]] * cls.N_SALT,
                            id=list(range(20, 20 + cls.N_SALT)))
        cls.system.part.add(type=[cls.types["Cl"]] * cls.N_SALT,
                            pos=np.random.rand(cls.N_SALT, 3) * cls.BOX_LENGTH,
                            q=[cls.charges["Cl"]] * cls.N_SALT)
        cls.system.part.add(type=[cls.types["HA"]] * cls.N_acid,
                            pos=np.random.rand(cls.N_acid, 3) * cls.BOX_LENGTH,
                            q=[cls.charges["HA"]] * cls.N_acid)

        cls.cph.add_reaction(
            gamma=10**(-cls.pKa),
            reactant_types=[cls.types["HA"]],
            product_types=[cls.types["A"], cls.types["H"]],
            default_charges={cls.types["HA"]: cls.charges["HA"],
                             cls.types["A"]: cls.charges["A"],
                             cls.types["H"]: cls.charges["H"]}
        )

        cls.widom.add_reaction(
            reactant_types=[],
            reactant_coefficients=[],
            product_types=[cls.types["Na"], cls.types["Cl"]],
            product_coefficients=[1, 1],
            default_charges={cls.types["Na"]: cls.charges["Na"],
                             cls.types["Cl"]: cls.charges["Cl"]}
        )
        cls.system.setup_type_map(type_list=list(cls.types.values()))

    def test_reaction_bookeeping(self):
        self.widom.calculate_particle_insertion_potential_energy(reaction_id=0)
        self.cph.reaction(reaction_steps=100)

        # Measure the chemical potential
        for _ in range(50):
            self.widom.calculate_particle_insertion_potential_energy(
                reaction_id=0)
            total_charge = sum(self.system.part.all().q)
            self.assertEqual(total_charge, 0.)


if __name__ == "__main__":
    ut.main()
