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
"""
This sample illustrates how particles of interest can be accessed via slicing.
"""
import espressomd
from espressomd import thermostat
import numpy as np

print("""
=======================================================
=                   slice_input.py                    =
=======================================================

Program Information:""")
print(espressomd.features())

# System parameters
#############################################################

box_l = 10.0

# Integration parameters
#############################################################

system = espressomd.System(box_l=[box_l] * 3)
system.set_random_state_PRNG()
#system.seed = system.cell_system.get_state()['n_nodes'] * [1234]
np.random.seed(seed=system.seed)

system.time_step = 0.01
system.cell_system.skin = 0.4

system.cell_system.max_num_cells = 2744


#############################################################
#  Setup System                                             #
#############################################################


# Particle setup
#############################################################

n_part = 10

id_list = np.arange(n_part)
pos_list = np.random.random((n_part, 3)) * system.box_l
type_list = np.ones(n_part, dtype=int)

system.part.add(id=id_list, pos=pos_list, type=type_list)


print("TYPE\n%s" % system.part[:].type)
system.part[0:2].type = [3, 3]
print("TYPE_NEW\n%s" % system.part[:].type)

print("POS\n%s" % system.part[:].pos)
system.part[:5].pos = [[1, 1, 1], [2, 2, 2], [3, 3, 3], [4, 4, 4], [5, 5, 5]]
print("POS_NEW\n%s" % system.part[:].pos)

print("V\n%s" % system.part[:].v)
system.part[:2].v = [[1, 2, 3], [2, 3, 4]]
print("V_NEW\n%s" % system.part[:].v)

print("F\n%s" % system.part[:].f)
system.part[:2].f = [[3, 4, 5], [4, 5, 6]]
print("F_NEW\n%s" % system.part[:].f)

if espressomd.has_features(["MASS"]):
    print("MASS\n%s" % system.part[:].mass)
    system.part[:2].mass = [2, 3]
    print("MASS_NEW\n%s" % system.part[:].mass)

if espressomd.has_features(["ELECTROSTATICS"]):
    print("Q\n%s" % system.part[:].q)
    system.part[::2].q = np.ones(n_part // 2)
    system.part[1::2].q = -np.ones(n_part // 2)
    print("Q_NEW\n%s" % system.part[:].q)
