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
"""
Illustrate how particles of interest can be accessed via slicing.
"""
import espressomd
import numpy as np

print("""
=======================================================
=                   slice_input.py                    =
=======================================================
""")

# System parameters
#############################################################

box_l = 10.0

# Integration parameters
#############################################################

system = espressomd.System(box_l=[box_l] * 3)
np.random.seed(seed=42)

system.time_step = 0.01
system.cell_system.skin = 0.4


#############################################################
#  Setup System                                             #
#############################################################


# Particle setup
#############################################################

n_part = 10

id_list = np.arange(n_part)
pos_list = np.random.random((n_part, 3)) * system.box_l
type_list = np.ones(n_part, dtype=int)

partcls = system.part.add(id=id_list, pos=pos_list, type=type_list)
p0p1 = system.part.by_ids([0, 1])

print("TYPE\n%s" % partcls.type)
p0p1.type = [3, 3]
print("TYPE_NEW\n%s" % partcls.type)

print("POS\n%s" % partcls.pos)
system.part.by_ids(range(5)).pos = [[1, 1, 1], [2, 2, 2], [
    3, 3, 3], [4, 4, 4], [5, 5, 5]]
print("POS_NEW\n%s" % partcls.pos)

print("V\n%s" % partcls.v)
p0p1.v = [[1, 2, 3], [2, 3, 4]]
print("V_NEW\n%s" % partcls.v)

print("F\n%s" % partcls.f)
p0p1.f = [[3, 4, 5], [4, 5, 6]]
print("F_NEW\n%s" % partcls.f)

if espressomd.has_features(["MASS"]):
    print("MASS\n%s" % partcls.mass)
    p0p1.mass = [2, 3]
    print("MASS_NEW\n%s" % partcls.mass)

if espressomd.has_features(["ELECTROSTATICS"]):
    print("Q\n%s" % partcls.q)
    system.part.by_ids(range(0, n_part, 2)).q = np.ones(n_part // 2)
    system.part.by_ids(range(1, n_part, 2)).q = -np.ones(n_part // 2)
    print("Q_NEW\n%s" % partcls.q)
