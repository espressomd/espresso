# Copyright (C) 2010-2019 The ESPResSo project
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


def WriteVTK(system, outFile):

    with open(outFile, "w") as fp:
        # header VTK
        fp.write("# vtk DataFile Version 2.0\n")
        fp.write("3D Unstructured Grid of Triangles\n")
        fp.write("ASCII\n")
        fp.write("DATASET UNSTRUCTURED_GRID\n")

        # points, get number from tables
        with open("tables/softPositions", "r") as fp2:
            numPoints = int(fp2.readline())
        fp.write(f"POINTS {numPoints} floats\n")

        # points, positions
        for p in system.part:
            fp.write(f"{' '.join(map(str, p.pos_folded))}\n")
