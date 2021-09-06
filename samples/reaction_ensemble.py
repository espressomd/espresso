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
"""
Guide for the reaction ensemble and the constant pH ensemble. The modeled
reaction is :math:`\\mathrm{AH} \\leftrightarrow \\mathrm{A}^- + \\mathrm{H}^+`.
"""
epilog = __doc__.split(":math:")[0] + "AH <-> A- + H+. " + \
    """You can choose in which ensemble you want to simulate via either
providing --reaction_ensemble or --constant_pH_ensemble as command line
argument to the script. Be aware that in the case of the reaction ensemble,
the dissociation constant gamma is not the thermodynamic reaction constant K,
but rather K * 1 mol/l and therefore carries a unit! In the case of the
constant pH method, gamma is the thermodynamic reaction constant!
"""
import numpy as np
import argparse

import espressomd
import espressomd.reaction_ensemble

parser = argparse.ArgumentParser(epilog=epilog)
group = parser.add_mutually_exclusive_group()
group.add_argument('--reaction_ensemble', action='store_const', dest='mode',
                   const='reaction_ensemble')
group.add_argument('--constant_pH_ensemble', action='store_const', dest='mode',
                   const='constant_pH_ensemble')
args = parser.parse_args()


# System parameters
#############################################################
box_l = 35

# Integration parameters
#############################################################
system = espressomd.System(box_l=[box_l] * 3)
np.random.seed(seed=42)

system.time_step = 0.02
system.cell_system.skin = 0.4

# Particle setup
#############################################################
types = {
    "HA": 0,
    "A-": 1,
    "H+": 2,
}
charge_dict = {
    types["HA"]: 0,
    types["A-"]: -1,
    types["H+"]: +1,
}

N0 = 50  # number of titratable units
K_diss = 0.0088

for i in range(N0):
    system.part.add(pos=np.random.random(3) * system.box_l, type=1)
for i in range(N0, 2 * N0):
    system.part.add(pos=np.random.random(3) * system.box_l, type=2)

RE = None
if args.mode == "reaction_ensemble":
    RE = espressomd.reaction_ensemble.ReactionEnsemble(
        kT=1,
        exclusion_radius=1,
        seed=77)
    RE.add_reaction(gamma=K_diss,
                    reactant_types=[types["HA"]],
                    reactant_coefficients=[1],
                    product_types=[types["A-"], types["H+"]],
                    product_coefficients=[1, 1],
                    default_charges=charge_dict)
elif args.mode == "constant_pH_ensemble":
    RE = espressomd.reaction_ensemble.ConstantpHEnsemble(
        kT=1, exclusion_radius=1, seed=77)
    RE.constant_pH = 2
    RE.add_reaction(gamma=K_diss, reactant_types=[types["HA"]],
                    product_types=[types["A-"], types["H+"]],
                    default_charges=charge_dict)
else:
    raise RuntimeError(
        "Please provide either --reaction_ensemble or --constant_pH_ensemble as argument ")

print(RE.get_status())
system.setup_type_map(list(types.values()))


# Set the hidden particle type to the lowest possible number to speed
# up the simulation
RE.set_non_interacting_type(max(types.values()) + 1)

for i in range(10000):
    RE.reaction()
    if i % 100 == 0:
        print("HA", system.number_of_particles(type=types["HA"]),
              "A-", system.number_of_particles(type=types["A-"]),
              "H+", system.number_of_particles(type=types["H+"]))

print("reaction 0 has acceptance rate: ", RE.get_acceptance_rate_reaction(0))
print("reaction 1 has acceptance rate: ", RE.get_acceptance_rate_reaction(1))
