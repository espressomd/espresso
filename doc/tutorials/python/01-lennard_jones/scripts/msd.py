#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
from espressomd import integrate
import cPickle as pickle
import os 
import numpy as np

data_file="data/sim_info.pickle"
part_file="data/part.pickle"

system = espressomd.System()
pickle.load(open(data_file,"r"))

print("""
ro variables:
cell_grid     {0.cell_grid}
cell_size     {0.cell_size} 
local_box_l   {0.local_box_l} 
max_cut       {0.max_cut}
max_part      {0.max_part}
max_range     {0.max_range} 
max_skin      {0.max_skin}
n_nodes       {0.n_nodes}
n_part        {0.n_part}
n_part_types  {0.n_part_types}
periodicity   {0.periodicity}
transfer_rate {0.transfer_rate}
""".format(system))

pickle.load(open(part_file,"r"))
print(system.part)
