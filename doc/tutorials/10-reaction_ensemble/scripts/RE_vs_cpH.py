#
# Copyright (C) 2013-2018 The ESPResSo project
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
from __future__ import print_function
from __future__ import division
import numpy as np

import espressomd
espressomd.assert_features(["ELECTROSTATICS", "LENNARD_JONES"])
from espressomd import code_info
from espressomd import analyze
from espressomd import integrate
from espressomd.interactions import *
from espressomd import reaction_ensemble


# System parameters
#
box_l = 60

# Integration parameters
#
system = espressomd.System(box_l=[box_l] * 3)
system.set_random_state_PRNG()
# system.seed = system.cell_system.get_state()['n_nodes'] * [1234]
np.random.seed(seed=12)

system.time_step = 0.02
system.cell_system.skin = 0.4
system.cell_system.max_num_cells = 2744


#
# Setup System                                             #
#
mode = "reaction_ensemble"
# mode="constant_pH_ensemble"

# Particle setup
#
# type 0 = HA
# type 1 = A-
# type 2 = H+

N0 = 50  # number of titratable units
K_diss = 0.000830

for i in range(N0):
    system.part.add(id=i, pos=np.random.random(3) * system.box_l, type=1)
for i in range(N0, 2 * N0):
    system.part.add(id=i, pos=np.random.random(3) * system.box_l, type=2)

RE = None
if(mode == "reaction_ensemble"):
    RE = reaction_ensemble.ReactionEnsemble(
        temperature=1,
        exclusion_radius=1,
     seed=2)
elif(mode == "constant_pH_ensemble"):
    RE = reaction_ensemble.ConstantpHEnsemble(
        temperature=1, exclusion_radius=1, seed=2)
    RE.constant_pH = 3
RE.add_reaction(gamma=K_diss, reactant_types=[0], reactant_coefficients=[1],
                product_types=[1, 2], product_coefficients=[1, 1],
                default_charges={0: 0, 1: -1, 2: +1})
print(RE.get_status())
system.setup_type_map([0, 1, 2])


alpha = []


for i in range(100):
    RE.reaction(100)
    print("HA", system.number_of_particles(type=0), "A-",
          system.number_of_particles(type=1), "H+",
          system.number_of_particles(type=2))
    alpha.append(system.number_of_particles(type=1) / N0)

alpha_av = np.mean(alpha)
alpha_err = np.std(alpha) / np.sqrt(len(alpha))

print("\n<alpha> = {} (err = {})".format(alpha_av, alpha_err))
