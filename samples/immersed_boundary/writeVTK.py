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
        fp.write("POINTS " + str(numPoints) + " floats\n")

        # points, positions
        for i in range(0, len(system.part)):
            fp.write(str(system.part[i].pos_folded[0]) + " " + str(
                system.part[i].pos_folded[1]) + " " + str(system.part[i].pos_folded[2]))
            fp.write("\n")
