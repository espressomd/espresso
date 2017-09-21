from __future__ import print_function
import espressomd
from espressomd.interactions import *

S = espressomd.System()

h = HarmonicBond(r_0=0, k=1)
f = FeneBond(k=1, d_r_max=1)
f2 = FeneBond(k=2, d_r_max=1)

print("\n**Defined three bond types:")
print(h)
print(f)
print(f2)

# note the order of the printed bond types
S.bonded_inter[0] = h
S.bonded_inter[2] = f
S.bonded_inter.add(f2)

print("\n**Added bond types to the system:")
for b in S.bonded_inter:
    print(b)

S.part.add(id=0, pos=(0., 0., 0.))
S.part.add(id=1, pos=(0, 0, 0))
S.part.add(id=2, pos=(0, 0, 0))
print("\n**Defined three particles")
print("Bonds for particle '0':")
print(S.part[0].bonds)

print("\n**Bonding particle 0 to particle 1 with bond refered by \"f2\'")
S.part[0].add_bond((f2, 1))
print("Bonds for particle '0':")
print(S.part[0].bonds)


print("\n**Bonding particle 0 to particle 2 with bond refered by index 0")
S.part[0].add_bond((0, 2))
print("Bonds for particle 0:")
print(S.part[0].bonds)

print("Bonds for particle 1:")
print(S.part[1].bonds)
print("Bonds for particle 2:")
print(S.part[2].bonds)


print("\n**creating holder for bonds on particle 0")
tmp = S.part[0].bonds
print("\n**deleted all bonds of 0")
S.part[0].delete_all_bonds()
print("Bonds for particle 0:")
print(S.part[0].bonds)

print("\n**setting particle 0 bonds from holder")
S.part[0].bonds = tmp
print("Bonds for particle 0:")
print(S.part[0].bonds)

print("\n**deleting bond refered by \"h\" to particle 2 :")
S.part[0].delete_bond((h, 2))
print("Bonds for particle '0':")
print(S.part[0].bonds)

print("\n**Bonding particle 0 to particle 1 with bondd refered by index 0 :")
print("**Bonding particle 0 to particle 2 with bondd refered by index 2 :")
S.part[0].bonds = ((0, 1), (2, 2))
print("Bonds for particle '0':")
print(S.part[0].bonds)

print("**Listing bonds held by all particles")
for p in S.part:
    print(p.bonds)
