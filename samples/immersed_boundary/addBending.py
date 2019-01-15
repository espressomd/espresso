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
from __future__ import print_function


def AddBending(system, kb):

    # currently only works for ONE SINGLE soft object

    # angles
    from espressomd.interactions import IBM_Tribend
    with open("tables/softAngles", "r") as fp:
        numAngles = int(fp.readline())
        print("Found {}".format(numAngles))
        # actual add
        for i in range(0, numAngles):
            line = str.split(fp.readline())
            id1 = int(line[0])
            id2 = int(line[1])
            id3 = int(line[2])
            id4 = int(line[3])
            tribend = IBM_Tribend(
                ind1=id1, ind2=id2, ind3=id3, ind4=id4, kb=kb, refShape="initial")
            system.bonded_inter.add(tribend)
            system.part[id1].add_bond((tribend, id2, id3, id4))
