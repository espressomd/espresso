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
import espressomd._system as es
import espressomd
from espressomd import code_info
from espressomd import analyze
from espressomd import integrate
from espressomd.interactions import *
from espressomd import reaction_ensemble
import numpy as np

import sys

print(code_info.features())
dev = "cpu"

# System parameters
#############################################################
box_l = 10

# Integration parameters
#############################################################
system = espressomd.System()
system.time_step = 0.02
system.skin = 0.4
system.max_num_cells = 2744


#############################################################
#  Setup System                                             #
#############################################################

# Interaction setup
#############################################################
system.box_l = [box_l, box_l, box_l]

# Particle setup
#############################################################
#type 0 = HA
#type 1 = A-
#type 2 = H+

N0 = 5 # number of titratable units
K_diss=0.0088

system.part.add(id=np.arange(N0) ,pos=np.random.random((N0,3)) * system.box_l, type=np.ones(N0)*1)
system.part.add(id=np.arange(N0,2*N0) ,pos=np.random.random((N0,3)) * system.box_l, type=2*np.ones(N0))


RE=reaction_ensemble.ReactionEnsemble(standard_pressure=0.00108, temperature=1, exclusion_radius=1)
RE.add(equilibrium_constant=K_diss,educt_types=[0],educt_coefficients=[1], product_types=[1,2], product_coefficients=[1,1])
RE.add(equilibrium_constant=1.0/K_diss, educt_types=[1,2], educt_coefficients=[1,1], product_types=[0],product_coefficients=[1])
RE.default_charges(dictionary={"0":0,"1":-1, "2":+1})
RE.print_status()

print("poss", system.part[:].pos)
while True:
	RE.reaction()
	_types=np.array(system.part[:].type)
	print("HA", len(_types[_types==0]), "A-", len(_types[_types==1]), "H+", len(_types[_types==2]))
	
