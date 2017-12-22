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
import numpy as np

import espressomd
from espressomd import code_info
from espressomd import analyze
from espressomd import integrate
from espressomd.interactions import *
from espressomd import reaction_ensemble
from espressomd import electrostatics

# System parameters
#############################################################
box_l = 10

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
    system.part.add(id=i, pos=np.random.random(3) * system.box_l, type=1, q=-1)
for i in range(N0, 2 * N0):
    system.part.add(id=i, pos=np.random.random(3) * system.box_l, type=2, q=1)

p3m = electrostatics.P3M(prefactor=0.5, accuracy=1e-2)
system.actors.add(p3m)
p3m_params = p3m.get_params()
for key in list(p3m_params.keys()):
    print("{} = {}".format(key, p3m_params[key]))
p3m.Tune(accuracy=1e3)

RE = reaction_ensemble.WidomInsertion(temperature=1, exclusion_radius=1.0)
RE.add(equilibrium_constant=K_diss, reactant_types=[], reactant_coefficients=[], product_types=[1, 2], product_coefficients=[1, 1])
RE.set_default_charges(dictionary={"1": -1, "2": +1})
print(RE.get_status())
system.setup_type_map([0, 1, 2])

for i in range(2000):
    RE.measure_excess_chemical_potential(1) #1 for insertion reaction
    for j in range(N0):
        RE.global_mc_move_for_one_particle_of_type(1)
        RE.global_mc_move_for_one_particle_of_type(2)
    if(i % 100 == 0):
        print("HA", system.number_of_particles(type=0), "A-",
              system.number_of_particles(type=1), "H+", system.number_of_particles(type=2))
              
print(RE.measure_excess_chemical_potential(1)) #1 for insertion reaction
