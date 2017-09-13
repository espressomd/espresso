from __future__ import print_function
import espressomd
from espressomd.interactions import *
from samples_common import open

system = espressomd.System()
f = FeneBond(k=1, d_r_max=1)
f2 = FeneBond(k=2, d_r_max=1.5)
h = HarmonicBond(r_0=0, k=1)

# Pickle data
###########################################################
try:
    import cPickle as pickle
except ImportError:
    import pickle

system.bonded_inter.add(f)
system.bonded_inter.add(f2)
system.bonded_inter.add(h)

output_filename = "bonded_inter_save"

with open(output_filename, "w") as bonded_ia_save:
    pickle.dump(system.bonded_inter, bonded_ia_save, -1)

print("The following bonding interactions were stored in file '{}':".format(
    output_filename))
for i in system.bonded_inter:
    print(i)
