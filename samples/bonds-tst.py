# Copyright (C) 2010-2018 The ESPResSo project
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
"""
This sample demonstrates how new bonds can be added between particles and how existing bond between particles can be deleted.
"""

from __future__ import print_function
import espressomd
from espressomd.interactions import *
from espressomd import assert_features

system = espressomd.System(box_l=[10.0, 10.0, 10.0])
system.seed = system.cell_system.get_state()['n_nodes'] * [1234]

h = HarmonicBond(r_0=0, k=1)
f = FeneBond(k=1, d_r_max=1)
f2 = FeneBond(k=2, d_r_max=1)

print("\n**Defined three bond types:")
print(h)
print(f)
print(f2)

# note the order of the printed bond types
system.bonded_inter[0] = h
system.bonded_inter[2] = f
system.bonded_inter.add(f2)

print("\n**Added bond types to the system:")
for b in system.bonded_inter:
    print(b)

system.part.add(id=0, pos=(0., 0., 0.))
system.part.add(id=1, pos=(0, 0, 0))
system.part.add(id=2, pos=(0, 0, 0))
print("\n**Defined three particles")
print("Bonds for particle '0':")
print(system.part[0].bonds)

print("\n**Bonding particle 0 to particle 1 with bond referred by \"f2\'")
system.part[0].add_bond((f2, 1))
print("Bonds for particle '0':")
print(system.part[0].bonds)


print("\n**Bonding particle 0 to particle 2 with bond referred by index 0")
system.part[0].add_bond((0, 2))
print("Bonds for particle 0:")
print(system.part[0].bonds)

print("Bonds for particle 1:")
print(system.part[1].bonds)
print("Bonds for particle 2:")
print(system.part[2].bonds)


print("\n**creating holder for bonds on particle 0")
tmp = system.part[0].bonds
print("\n**deleted all bonds of 0")
system.part[0].delete_all_bonds()
print("Bonds for particle 0:")
print(system.part[0].bonds)

print("\n**setting particle 0 bonds from holder")
system.part[0].bonds = tmp
print("Bonds for particle 0:")
print(system.part[0].bonds)

print("\n**deleting bond referred by \"h\" to particle 2 :")
system.part[0].delete_bond((h, 2))
print("Bonds for particle '0':")
print(system.part[0].bonds)

print("\n**Bonding particle 0 to particle 1 with bond referred by index 0 :")
print("**Bonding particle 0 to particle 2 with bond referred by index 2 :")
system.part[0].bonds = ((0, 1), (2, 2))
print("Bonds for particle '0':")
print(system.part[0].bonds)

print("**Listing bonds held by all particles")
for p in system.part:
    print(p.bonds)
