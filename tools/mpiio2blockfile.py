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

# Important: Currently, this script assumes, that all fields are dumped.
# To support less: Read the "fields" unsigned int from 1.head.

import os
import sys
import array

if len(sys.argv) < 2:
    print >>sys.stderr, "Prefix as second argument!"
    raise SystemExit(1)

fpref = sys.argv[1]
headf = fpref + ".head"
preff = fpref + ".pref"
idf = fpref + ".id"
typef = fpref + ".type"
posf = fpref + ".pos"
velf = fpref + ".vel"
bofff = fpref + ".boff"
bondf = fpref + ".bond"

# Determine nproc, etc. at time of writing
nproc = os.stat(preff).st_size / 4
ntotalpart = os.stat(idf).st_size / 4
ntotalbond = os.stat(bondf).st_size / 4

# Read header - fields, n_bonded_ia, bonded_ia_params[:].nums
fields = 0
nbia = 0
biaparams = array.array("i")
with open(headf) as f:
    afields = array.array("I")
    afields.read(f, 1)
    fields = afields[0]
    anbia = array.array("i")
    anbia.read(f, 1)
    nbia = anbia[0]
    biaparams.read(f, nbia)

# Read prefixes
pref = array.array("i")
with open(preff) as f:
    pref.read(f, nproc)

# Read ids
id = array.array("i")
with open(idf) as f:
    id.read(f, ntotalpart)

# Read types
type = array.array("i")
with open(typef) as f:
    type.read(f, ntotalpart)

# Read pos
pos = array.array("d")
with open(posf) as f:
    pos.read(f, 3 * ntotalpart)

# Read vel
vel = array.array("d")
with open(velf) as f:
    vel.read(f, 3 * ntotalpart)

# Read bonds
boff = array.array("i")
with open(bofff) as f:
    boff.read(f, ntotalpart + nproc)

bond = array.array("i")
with open(bondf) as f:
    bond.read(f, ntotalbond)


# Print particles in blockfile format
print("{particles {id type pos v}")
for i in xrange(ntotalpart):
    print("\t{%i %i %r %r %r %r %r %r}" \
        % (id[i], type[i], pos[3*i], pos[3*i+1], pos[3*i+2], \
           vel[3*i], vel[3*i+1], vel[3*i+2]))
print("}")

# Print bonds in blockfile format
print("{bonds")
addend = 0 # ntotal bonds of previous processors
for rank in xrange(nproc):
    # The start and end indices for the boff array are determined via
    # pref. However, there are (nlocalpart + 1) boff entries per proc.
    start = pref[rank] + rank
    end = rank + (pref[rank + 1] if rank < nproc - 1 else ntotalpart)
    for pid, i in enumerate(xrange(start, end)):
        print("\t{%i { " % id[pref[rank] + pid], end="")
        # The start and end indices for the bond array are determined
        # via boff. However, boff does only *locally* store prefixes,
        # i.e. they have to be globalized by adding the total number of
        # bonds on all ranks before this one.
        j = addend + boff[i]
        while j < addend + boff[i + 1]:
            bond_num = bond[j]
            j += 1
            npartners = biaparams[bond_num]
            print("{%i" % bond_num, end="")
            for _ in xrange(npartners):
                print(" %i" % bond[j], end="")
                j += 1
            print("} ", end="")
        print("} }")
    addend += boff[end]
print("}")
