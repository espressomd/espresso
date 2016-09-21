from __future__ import print_function
import espressomd
from espressomd.interactions import *

try:
    import cPickle as pickle
except ImportError:
    import pickle

# Import system properties
#############################################################
system = espressomd.System()

# demo: all previously defined bonded-interactions will be 
# overwritten intentionally after pickle.load such that the
# original state from pickle.dump will be restored
f_new = FeneBond(k=3, d_r_max=2)
h_new = HarmonicBond(r_0=5, k=4)
system.bonded_inter.add(f_new)
system.bonded_inter.add(h_new)

print("Initial bonded interactions:")
for i in system.bonded_inter:
    print(i)

with open("bonded_inter_save", "r") as bonded_ia_save:
    pickle.load(bonded_ia_save)

print("\n\nBonded interactions after loading:")
for i in system.bonded_inter:
    print(i)
