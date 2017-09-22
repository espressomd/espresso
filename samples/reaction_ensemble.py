#
# Copyright (C) 2013,2014 The ESPResSo project
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
import espressomd
from espressomd import code_info
from espressomd import analyze
from espressomd import integrate
from espressomd.interactions import *
from espressomd import reaction_ensemble
from espressomd import grand_canonical
import numpy as np

import sys

if('REACTION_ENSEMBLE' not in espressomd.code_info.features()):
    print("REACTION_ENSEMBLE not compiled in.")
    sys.exit()
dev = "cpu"

# System parameters
#############################################################
box_l = 35

# Integration parameters
#############################################################
system = espressomd.System()
system.time_step = 0.02
system.cell_system.skin = 0.4
system.cell_system.max_num_cells = 2744


#############################################################
#  Setup System                                             #
#############################################################

# Interaction setup
#############################################################
system.box_l = [box_l, box_l, box_l]

# Particle setup
#############################################################
# type 0 = HA
# type 1 = A-
# type 2 = H+

N0 = 50  # number of titratable units
K_diss = 0.0088

for i in range(N0):
    system.part.add(id=i, pos=np.random.random(3) * system.box_l, type=1)
for i in range(N0, 2 * N0):
    system.part.add(id=i, pos=np.random.random(3) * system.box_l, type=2)


RE = reaction_ensemble.reaction_ensemble(
    standard_pressure=0.00108, temperature=1, exclusion_radius=1)
RE.add(equilibrium_constant=K_diss, reactant_types=[0], reactant_coefficients=[
       1], product_types=[1, 2], product_coefficients=[1, 1])
RE.set_default_charges(dictionary={"0": 0, "1": -1, "2": +1})
print(RE.get_status())
grand_canonical.setup([0, 1, 2])

# RE.pH=2

for i in range(10000):
    RE.reaction()
#	RE.reaction_constant_pH()
    if(i % 100 == 0):
        print("HA", grand_canonical.number_of_particles(current_type=0), "A-",
              grand_canonical.number_of_particles(current_type=1), "H+", grand_canonical.number_of_particles(current_type=2))
