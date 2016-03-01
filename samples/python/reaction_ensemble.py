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
from espressomd import thermostat
from espressomd import code_info
from espressomd import analyze
from espressomd import integrate
from espressomd import reaction_ensemble
import numpy as np

import sys

print("""
=======================================================
=                   slice_input.py                    =
=======================================================

Program Information:""")
print(code_info.features())

dev = "cpu"

# System parameters
#############################################################

box_l = 10.0

# Integration parameters
#############################################################

system = espressomd.System()
system.time_step = 0.01
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

n_part = 10

id_list = np.arange(n_part)
pos_list = np.random.random((n_part,3)) * system.box_l
type_list = [0,0,0,0,0,1,1,1,1,1]

system.part.add(id=id_list ,pos=pos_list, type=type_list)


K_diss=0.0088

RE=reaction_ensemble.ReactionEnsemble(standard_pressure=0.00108, temperature=1, exclusion_radius=1)
RE.add(equilibrium_constant=K_diss,educt_types=[0,1],educt_coefficients=[1,1], product_types=[2], product_coefficients=[1])
RE.add(equilibrium_constant=1.0/K_diss,product_types=[0,1],product_coefficients=[1,1], educt_types=[2], educt_coefficients=[1])
RE.default_charges(dictionary={"0":1,"1":-1, "2":0})

RE.fix_polymer_monomers=True

for i in range(100000) :
	_types=np.array(system.part[:].type)
	RE.reaction()
	_types=np.array(system.part[:].type)
#	print(len(_types[_types==0]))

