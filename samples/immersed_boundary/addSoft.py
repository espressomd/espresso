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


def AddSoft(system, comX, comY, comZ, k1, k2):

    # currently only works for ONE SINGLE soft object

    # open file and add nodes
    with open("tables/softPositions", "r") as fp:
        numPoints = int(fp.readline())
        print("Found {} nodes".format(numPoints))

        # actual add
        for i in range(0, numPoints):
            line = str.split(fp.readline())
            X = float(line[0]) + comX
            Y = float(line[1]) + comY
            Z = float(line[2]) + comZ
#            print X, Y, Z
            system.part.add(id=i, pos=[X, Y, Z], virtual=1)

    # triangles
    from espressomd.interactions import IBM_Triel
    with open("tables/softTriangles", "r") as fp:
        numTri = int(fp.readline())
        print("Found {} triangles".format(numTri))
        # actual add
        for i in range(0, numTri):
            line = str.split(fp.readline())
            id1 = int(line[0])
            id2 = int(line[1])
            id3 = int(line[2])
            tri = IBM_Triel(ind1=id1, ind2=id2, ind3=id3,
                            elasticLaw="Skalak", k1=k1, k2=k2, maxDist=5)
            system.bonded_inter.add(tri)
            system.part[id1].add_bond((tri, id2, id3))
