from __future__ import print_function
import espressomd
from espressomd.interactions import *

S = espressomd.System()
h = HarmonicBond(r_0=0, k=1)
f = FeneBond(k=1, d_r_max=1)
f2 = FeneBond(k=2, d_r_max=1)
print(f2)


S.bonded_inter[0] = h
S.bonded_inter[2] = f
S.bonded_inter.add(f2)

for b in S.bonded_inter:
    print(b)

S.part.add(id=0, pos=(0., 0., 0.))
S.part.add(id=1, pos=(0, 0, 0))
S.part.add(id=2, pos=(0, 0, 0))


print(S.part[0].bonds)
S.part[0].add_bond((f2, 1))
S.part[0].add_bond((0, 2))
print(S.part[0].bonds)
tmp = S.part[0].bonds
S.part[0].delete_all_bonds()
print(S.part[0].bonds)
S.part[0].bonds = tmp
print(S.part[0].bonds)
S.part[0].delete_bond((h, 2))
print(S.part[0].bonds)
S.part[0].bonds = ((0, 1), (2, 2))
print(S.part[0].bonds)
