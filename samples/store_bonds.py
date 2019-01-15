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
This sample illustrates how bond information can be stored.
"""

from __future__ import print_function
import espressomd
from espressomd import interactions

system = espressomd.System(box_l=[10.0, 10.0, 10.0])
f = interactions.FeneBond(k=1, d_r_max=1)
f2 = interactions.FeneBond(k=2, d_r_max=1.5)
h = interactions.HarmonicBond(r_0=0, k=1)

# Pickle data
###########################################################
try:
    import cPickle as pickle
except ImportError:
    import pickle

system.bonded_inter.add(f)
system.bonded_inter.add(f2)
system.bonded_inter.add(h)

output_filename = "bonded_inter_save.pkl"

with open(output_filename, "wb") as bonded_ia_save:
    pickle.dump(system.bonded_inter, bonded_ia_save, -1)

print("The following bonding interactions were stored in file '{}':".format(
    output_filename))
for i in system.bonded_inter:
    print(i)
