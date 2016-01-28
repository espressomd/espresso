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
import numpy as np

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

for i in range(10):
    system.part.add(id=i, pos=[0,0,0], type=0)

print(system.part[0:2].type)
system.part[0:2].type=[3,3]
print(system.part[:].type)

print(system.part[0:3].pos)
system.part[::2].pos=[[1,1,1],[2,2,2],[3,3,3],[4,4,4],[5,5,5]]
print(system.part[:].pos)

exit()

n_part = 1000
id_list=np.arange(n_part)
pos_list=np.random.random((n_part,3)) * system.box_l

print(id_list)
print(pos_list)

system.part.add(id=id_list,pos=pos_list)

# Assingn charge to particles

system.part[0:n_part:2].q = -1.0
system.part[1:n_part:2].q = 1.0



