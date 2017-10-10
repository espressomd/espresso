import numpy as np
from oif_utils import *
from espressomd.interactions import OifLocalForces
from espressomd.interactions import OifGlobalForces
from espressomd.interactions import OifOutDirection



class FixedPoint:
    """ Represents mesh points, not connected to any ESPResSo particle"""
    def __init__(self, pos, id):
        self.x = pos[0]
        self.y = pos[1]
        self.z = pos[2]
        self.id = id

    def GetPos(self):
        return [self.x, self.y, self.z]

    def GetID(self):
        return self.id

    # def SetPos(self, pos):
    #     self.x = pos[0]
    #     self.y = pos[1]
    #     self.z = pos[2]


class PartPoint:
    """ Represents mesh points, connected to ESPResSo particle"""
    def __init__(self, part, id, partId):  # part is physical ESPResSo particle corresponding to that particular point
        self.part = part
        self.partId = partId  # because in adding bonds to the particles in OifCell one needs to know the global id of the particle.    
        self.id = id

    def GetPos(self):
        return self.part.pos

    def GetVel(self):
        return self.part.v

    def GetMass(self):
        return self.part.mass

    def GetType(self):
        return self.part.type

    def GetForce(self):
        return self.part.f

    def SetPos(self,pos):
        self.part.pos = pos

    def SetVel(self, vel):
        self.part.v = vel

    def SetForce(self, force):
        self.part.ext_force = force

    def KillMotion(self):
        self.part.fix = [1, 1, 1]
        
    def UnkillMotion(self):
        self.part.unfix()


class Edge:
    """ Represents edges in a mesh"""

    def __init__(self, A, B):
        self.A = A
        self.B = B

    def Length(self):
        return Distance(self.A.GetPos(), self.B.GetPos())


class Triangle:
    """ Represents triangles in a mesh"""

    def __init__(self, A, B, C):
        self.A = A
        self.B = B
        self.C = C

    def Area(self):
        area = AreaTriangle(self.A.GetPos(), self.B.GetPos(), self.C.GetPos())
        return area


class Angle:
    """ Represents triangles in a mesh"""

    def __init__(self, A, B, C, D):
        self.A = A
        self.B = B
        self.C = C
        self.D = D

    def Size(self):
        angleSize = AngleBtwTriangles(self.A.GetPos(), self.B.GetPos(), self.C.GetPos(), self.D.GetPos())
        return angleSize

class ThreeNeighbors:
    """ Represents three best spatially distributed neighbors of a point in a mesh"""

    def __init__(self, A, B, C):
        self.A = A
        self.B = B
        self.C = C

    def OuterNormal(self):
        outerNormal = GetNTriangle(self.A.GetPos(), self.B.GetPos(), self.C.GetPos()) #check this!!!!
        return outerNormal

class Mesh:
    """ Represents a triangular mesh """

    def __init__(self, nodesFile=None, trianglesFile=None, system=None, stretch=(1.0, 1.0, 1.0), partType=-1, partMass=1.0, normal=False, checkOrientation=True):
        if system is None:
            print "Mesh: No system provided. Quitting."
            quit()
        self.system = system
        self.normal = normal
        self.nodesFile = nodesFile
        self.trianglesFile = trianglesFile

        self.points = []
        self.edges = []
        self.triangles = []
        self.angles = []
        self.neighbors = []
        self.idsExtremalPoints = [0, 0, 0, 0, 0, 0, 0]
        
        if (nodesFile is None) or (trianglesFile is None):
            print "OifMesh warning: one of nodesFile or trianglesFile was not given. May raise concern when Mesh object was being created."
        else:
            # reading the mesh point positions from file
            infile = open(nodesFile, "r")
            nodesCoord = infile.read().split("\n")
            infile.close()
            # removes a blank line at the end of the file if there is any:
            nodesCoord = filter(None, nodesCoord) # here we have list of lines with triplets of strings
            for line in nodesCoord:  # extracts coordinates from the string line
                line = line.split()
                coords = [stretch[0]*float(line[0]), stretch[1]*float(line[1]), stretch[2]*float(line[2])]
                tmpFixedPoint = FixedPoint(coords, len(self.points))
                self.points.append(tmpFixedPoint)
            #print "Coordinates of the mesh points"
            #for tmpFixedPoint in self.points:
            #    print tmpFixedPoint.GetPos()

            # searching for extremal points IDs
            xMin = largeNumber
            xMax = -largeNumber
            yMin = largeNumber
            yMax = -largeNumber
            zMin = largeNumber
            zMax = -largeNumber
            for tmpFixedPoint in self.points:
                coords = tmpFixedPoint.GetPos()
                if coords[0] < xMin:
                    xMin = coords[0]
                    self.idsExtremalPoints[0] = tmpFixedPoint.GetID()
                if coords[0] > xMax:
                    xMax = coords[0]
                    self.idsExtremalPoints[1] = tmpFixedPoint.GetID()
                if coords[1] < yMin:
                    yMin = coords[1]
                    self.idsExtremalPoints[2] = tmpFixedPoint.GetID()
                if coords[1] > yMax:
                    yMax = coords[1]
                    self.idsExtremalPoints[3] = tmpFixedPoint.GetID()
                if coords[2] < zMin:
                    zMin = coords[2]
                    self.idsExtremalPoints[4] = tmpFixedPoint.GetID()
                if coords[2] > zMax:
                    zMax = coords[2]
                    self.idsExtremalPoints[5] = tmpFixedPoint.GetID()

            # reading the triangle incidences from file
            infile = open(trianglesFile, "r")
            trianglesIncid = infile.read().split("\n")
            infile.close()
            # removes a blank line at the end of the file if there is any:
            trianglesIncid = filter(None, trianglesIncid)
            for line in trianglesIncid:  # extracts incidences from the string line
                line = line.split()
                incid = [int(line[0]), int(line[1]), int(line[2])]
                tmpTriangle = Triangle(self.points[incid[0]], self.points[incid[1]], self.points[incid[2]])
                self.triangles.append(tmpTriangle)

            if checkOrientation is True:
                # check whether all triangles in file had the same orientation; if not, correct the orientation
                print "OifMesh: Checking orientation of triangles (and repairing if needed). " \
                      "For large meshes this may take a while. " \
                      "If you are certain your mesh is correct, this can be skipped using the OifCellType option checkOrientation = False."
                self.CheckOrientation()

            #print "Triangle incidences"
            #for tmpTriangle in self.triangles:
            #    print tmpTriangle.A.GetPos(), tmpTriangle.B.GetPos(), tmpTriangle.C.GetPos()

            # creating list of edge incidences from triangle incidences
            # using temporary list of edge incidences
            tmpEdgeIncidences = []
            for triangle in self.triangles:
                pa = triangle.A.id
                pb = triangle.B.id
                pc = triangle.C.id
                if (([pa, pb] not in tmpEdgeIncidences) and ([pb, pa] not in tmpEdgeIncidences)):
                    tmpEdgeIncidences.append([pa, pb])
                if (([pb, pc] not in tmpEdgeIncidences) and ([pc, pb] not in tmpEdgeIncidences)):
                    tmpEdgeIncidences.append([pb, pc])
                if (([pa, pc] not in tmpEdgeIncidences) and ([pc, pa] not in tmpEdgeIncidences)):
                    tmpEdgeIncidences.append([pa, pc])
            for tmpIncid in tmpEdgeIncidences:
                tmpEdge = Edge(self.points[tmpIncid[0]], self.points[tmpIncid[1]])
                self.edges.append(tmpEdge)

            #print "Edges:"
            #for tmpEdge in self.edges:
            #    print tmpEdge.A.id, tmpEdge.B.id

            ## creating list angles (former bending incidences) from triangle incidences
            for edge in self.edges:
                pa = edge.A.id
                pb = edge.B.id
                detected = 0
                # detected = number of detected triangles with current edge common
                # Algorithm is as follows: we run over all triangles and check whether two vertices are those from current edge. If we find such triangle, we put the ID of the third vertex to pc and we check if the orientation pa, pb, pc is the same as was in the triangle list (meaning, that we found one of the following three triples in the triangle list: pa, pb, pc or pb, pc, pa or pc, pa, pb). If we have the same orientation, we set orient = 1, otherwise orient = -1.
                # Then we go further looking for the second triangle. The second triangle should have the opposite orientation.
                # The normal of the first triangle will be P1P2 x P1P3, of the second triangle will be P2P4 x P2P3
                orient = 0
                for triangle in self.triangles:
                    # Run over all triangles and determine the two triangles with the common current edge
                    if ((pa == triangle.A.id) and (pb == triangle.B.id)):
                        if (detected == 0):
                            # if no triangle with such edge was detected before
                            pc = triangle.C.id
                            detected = 1
                            orient = 1
                        else:
                            # if this is the second triangle with this edge, then also quit the for-loop over triangles 
                            pd = triangle.C.id
                            break
                    if ((pa == triangle.B.id) and (pb == triangle.C.id)):
                        if (detected == 0):
                            pc = triangle.A.id
                            detected = 1
                            orient = 1
                        else:
                            pd = triangle.A.id
                            break
                    if ((pa == triangle.C.id) and (pb == triangle.A.id)):
                        if (detected == 0):
                            pc = triangle.B.id
                            detected = 1
                            orient = 1
                        else:
                            pd = triangle.B.id
                            break
                    if ((pa == triangle.B.id) and (pb == triangle.A.id)):
                        if (detected == 0):
                            pc = triangle.C.id
                            detected = 1
                            orient = -1
                        else:
                            pd = triangle.C.id
                            break
                    if ((pa == triangle.C.id) and (pb == triangle.B.id)):
                        if (detected == 0):
                            pc = triangle.A.id
                            detected = 1
                            orient = -1
                        else:
                            pd = triangle.A.id
                            break
                    if ((pa == triangle.A.id) and (pb == triangle.C.id)):
                        if (detected == 0):
                            pc = triangle.B.id
                            detected = 1
                            orient = -1
                        else:
                            pd = triangle.B.id
                            break
                if (orient == 1):
                    tmp = pd
                    pd = pc
                    pc = tmp
                tmpAngle = Angle(self.points[pc], self.points[pa], self.points[pb], self.points[pd])
                self.angles.append(tmpAngle)

            #print "Angles:"
            #for tmpAngle in self.angles:
            #    print tmpAngle.A.id, tmpAngle.B.id, tmpAngle.C.id, tmpAngle.D.id

            # creating list of three neighbors for membrane collision
            if normal is True:
                for point in self.points:
                    tmpNeighbors = []
                    # cycle through edges and select those that contain point
                    for edge in self.edges:
                        # take an edge and copy the nodes of the edge to pa, pb
                        if edge.A.id == point.id:
                            tmpNeighbors.append(edge.B)
                        if edge.B.id == point.id:
                            tmpNeighbors.append(edge.A)
                    # create vectors to all neighbors and normalize them
                    tmpVectorsToNeighbors = []
                    pCoords = np.array(point.GetPos())
                    for neighbor in tmpNeighbors:
                        tmpVector = neighbor.GetPos() - pCoords
                        tmpLength = Norm(tmpVector)
                        if tmpLength < smallEpsilon:
                            print "Mesh: Degenerate edge. Quitting"
                            quit()
                        tmpVector /= tmpLength
                        tmpVectorsToNeighbors.append(tmpVector)
                    # check all triplets of neighbors and select the one that is best spatially distributed
                    # by adding the corresponding three normalized vectors
                    # and selecting the one with smallest resultant vector
                    nNeighbors = len(tmpNeighbors)
                    minLength = largeNumber
                    bestNeighbors = [tmpNeighbors[0], tmpNeighbors[1], tmpNeighbors[2]]
                    for i in range(0,nNeighbors):
                        for j in range(i+1,nNeighbors):
                            for k in range(j+1,nNeighbors):
                                tmpResultVector = tmpVectorsToNeighbors[i] + tmpVectorsToNeighbors[j] + tmpVectorsToNeighbors[k]
                                tmpResultVectorLength = Norm(tmpResultVector)
                                if tmpResultVectorLength < minLength:
                                    minLength = tmpResultVectorLength
                                    bestNeighbors = [tmpNeighbors[i], tmpNeighbors[j], tmpNeighbors[k]]
                    # find one triangle that contains this point and compute its normal vector
                    for triangle in self.triangles:
                        if (triangle.A.id == point.id or triangle.B.id == point.id or triangle.C.id == point.id):
                            tmpNormalTriangle = GetNTriangle(triangle.A.GetPos(), triangle.B.GetPos(), triangle.C.GetPos())
                            break
                    # properly orient selected neighbors and save them to the list of neighbors
                    tmpNormalNeighbors = GetNTriangle(bestNeighbors[0].GetPos(), bestNeighbors[1].GetPos(), bestNeighbors[2].GetPos())
                    tmpLengthNormalTriangle = Norm(tmpNormalTriangle)
                    tmpLengthNormalNeighbors = Norm(tmpNormalNeighbors)
                    tmpProduct = (tmpNormalTriangle[0] * tmpNormalNeighbors[0] + tmpNormalTriangle[1] * tmpNormalNeighbors[1] + tmpNormalTriangle[2] * tmpNormalNeighbors[2]) / (tmpLengthNormalTriangle * tmpLengthNormalNeighbors)
                    tmpAngle = np.arccos(tmpProduct)
                    if (tmpAngle > np.pi/2.0):
                        selectedNeighbors = ThreeNeighbors(bestNeighbors[0], bestNeighbors[1], bestNeighbors[2])
                    else:
                        selectedNeighbors = ThreeNeighbors(bestNeighbors[0], bestNeighbors[2], bestNeighbors[1])
                    self.neighbors.append(selectedNeighbors)
            else:
                for point in self.points:
                    selectedNeighbors = ThreeNeighbors(point, point, point)
                    self.neighbors.append(selectedNeighbors)

            #print "Neighbors:"
            #for tmpNeighbor in self.neighbors:
            #    print tmpNeighbor.A.id, tmpNeighbor.B.id, tmpNeighbor.C.id

    def Copy(self, origin=None, partType=-1, partMass=1.0, rotate=None):
        mesh = Mesh(system=self.system)
        mesh.idsExtremalPoints = self.idsExtremalPoints

        if rotate is not None:
            # variables for rotation
            ca = np.cos(rotate[0])
            sa = np.sin(rotate[0])
            cb = np.cos(rotate[1])
            sb = np.sin(rotate[1])
            cc = np.cos(rotate[2])
            sc = np.sin(rotate[2])
            rotation = np.array([[cb * cc, sa * sb * cc - ca * sc, sc * sa + cc * sb * ca],
                                [cb * sc, ca * cc + sa * sb * sc, sc * sb * ca - cc * sa],
                                [-sb, cb * sa, ca * cb]])

        for point in self.points:
            # PartPoints are created
            tmpPos = point.GetPos()
            tmpRotatePos = np.array(point.GetPos())
            # rotation of nodes
            if rotate is not None:
                tmpPos = rotation.dot(tmpRotatePos)
                tmpPos = [DiscardEpsilon(tmpPos[0]), DiscardEpsilon(tmpPos[1]), DiscardEpsilon(tmpPos[2])]
            if origin is not None:
                tmpPos[0] += origin[0]
                tmpPos[1] += origin[1]
                tmpPos[2] += origin[2]
            newPartId = len(self.system.part)  # to remember the global id of the ESPResSo particle
            self.system.part.add(pos=tmpPos, type=partType, mass=partMass, mol_id=partType)
            newPart = self.system.part[newPartId]
            newPartPoint = PartPoint(newPart, len(mesh.points), newPartId)
            mesh.points.append(newPartPoint)
        for edge in self.edges:
            newEdge = Edge(mesh.points[edge.A.id], mesh.points[edge.B.id])
            mesh.edges.append(newEdge)
        for triangle in self.triangles:
            newTriangle = Triangle(mesh.points[triangle.A.id], mesh.points[triangle.B.id], mesh.points[triangle.C.id])
            mesh.triangles.append(newTriangle)
        for angle in self.angles:
            newAngle = Angle(mesh.points[angle.A.id], mesh.points[angle.B.id], mesh.points[angle.C.id], mesh.points[angle.D.id])
            mesh.angles.append(newAngle)
        for neighbors in self.neighbors:
            newNeighbors = ThreeNeighbors(mesh.points[neighbors.A.id], mesh.points[neighbors.B.id], mesh.points[neighbors.C.id])
            mesh.neighbors.append(newNeighbors)
        # skopirovat aj normal true/false
        # skopiruje obsah....
        return mesh

    def CheckOrientation(self):
        tmpTriangleList = []
        tmpTriangleListOK = []
        for triangle in self.triangles:
            tmpTriangleList.append(triangle)

        # move the first triangle to the checked and corrected list
        tmpTriangleListOK.append(tmpTriangleList[0])
        tmpTriangleList.pop(0)

        #print strftime("%H:%M:%S", gmtime())
        while len(tmpTriangleList) is not 0:
            i = 0
            while i < len(tmpTriangleList):
                tmpTriangle = tmpTriangleList[i]
                for correctTriangle in tmpTriangleListOK:
                    # check if triangles have a common edge, if so, check orientation
                    areNeighbors = True
                    if tmpTriangle.A.id == correctTriangle.A.id:
                        if tmpTriangle.B.id == correctTriangle.B.id:
                            tOK = False  # this is situation 123 and 124
                            correctedTriangle = Triangle(tmpTriangle.A, tmpTriangle.C, tmpTriangle.B)
                        else:
                            if tmpTriangle.B.id == correctTriangle.C.id:
                                tOK = True  # this is situation 123 and 142
                            else:
                                if tmpTriangle.C.id == correctTriangle.B.id:
                                    tOK = True  # this is situation 123 and 134
                                else:
                                    if tmpTriangle.C.id == correctTriangle.C.id:
                                        tOK = False  # this is situation 123 and 143
                                        correctedTriangle = Triangle(tmpTriangle.A, tmpTriangle.C, tmpTriangle.B)
                                    else:
                                        areNeighbors = False
                    else:
                        if tmpTriangle.A.id == correctTriangle.B.id:
                            if tmpTriangle.B.id == correctTriangle.C.id:
                                tOK = False  # this is situation 123 and 412
                                correctedTriangle = Triangle(tmpTriangle.A, tmpTriangle.C, tmpTriangle.B)
                            else:
                                if tmpTriangle.B.id == correctTriangle.A.id:
                                    tOK = True  # this is situation 123 and 214
                                else:
                                    if tmpTriangle.C.id == correctTriangle.C.id:
                                        tOK = True  # this is situation 123 and 413
                                    else:
                                        if tmpTriangle.C.id == correctTriangle.A.id:
                                            tOK = False  # this is situation 123 and 314
                                            correctedTriangle = Triangle(tmpTriangle.A, tmpTriangle.C, tmpTriangle.B)
                                        else:
                                            areNeighbors = False
                        else:
                            if tmpTriangle.A.id == correctTriangle.C.id:
                                if tmpTriangle.B.id == correctTriangle.A.id:
                                    tOK = False  # this is situation 123 and 241
                                    correctedTriangle = Triangle(tmpTriangle.A, tmpTriangle.C, tmpTriangle.B)
                                else:
                                    if tmpTriangle.B.id == correctTriangle.B.id:
                                        tOK = True  # this is situation 123 and 421
                                    else:
                                        if tmpTriangle.C.id == correctTriangle.A.id:
                                            tOK = True  # this is situation 123 and 341
                                        else:
                                            if tmpTriangle.C.id == correctTriangle.B.id:
                                                tOK = False  # this is situation 123 and 431
                                                correctedTriangle = Triangle(tmpTriangle.A, tmpTriangle.C,
                                                                             tmpTriangle.B)
                                            else:
                                                areNeighbors = False
                            else:
                                if tmpTriangle.B.id == correctTriangle.A.id:
                                    if tmpTriangle.C.id == correctTriangle.B.id:
                                        tOK = False  # this is situation 123 and 234
                                        correctedTriangle = Triangle(tmpTriangle.A, tmpTriangle.C, tmpTriangle.B)
                                    else:
                                        if tmpTriangle.C.id == correctTriangle.C.id:
                                            tOK = True  # this is situation 123 and 243
                                        else:
                                            areNeighbors = False
                                else:
                                    if tmpTriangle.B.id == correctTriangle.B.id:
                                        if tmpTriangle.C.id == correctTriangle.C.id:
                                            tOK = False  # this is situation 123 and 423
                                            correctedTriangle = Triangle(tmpTriangle.A, tmpTriangle.C,
                                                                         tmpTriangle.B)
                                        else:
                                            if tmpTriangle.C.id == correctTriangle.A.id:
                                                tOK = True  # this is situation 123 and 324
                                            else:
                                                areNeighbors = False
                                    else:
                                        if tmpTriangle.B.id == correctTriangle.C.id:
                                            if tmpTriangle.C.id == correctTriangle.A.id:
                                                tOK = False  # this is situation 123 and 342
                                                correctedTriangle = Triangle(tmpTriangle.A, tmpTriangle.C,
                                                                             tmpTriangle.B)
                                            else:
                                                if tmpTriangle.C.id == correctTriangle.B.id:
                                                    tOK = True  # this is situation 123 and 432
                                                else:
                                                    areNeighbors = False
                                        else:
                                            areNeighbors = False
                    if areNeighbors:
                        # move the tmpTriangle to the checked and corrected list
                        if tOK:
                            tmpTriangleListOK.append(tmpTriangle)
                        else:
                            tmpTriangleListOK.append(correctedTriangle)
                            print "OifMesh: Correcting orientation of triangle."
                        tmpTriangleList.pop(i)
                        break
                i += 1
        # replace triangles with checked triangles
        i = 0
        for tmpTriangle in tmpTriangleListOK:
            self.triangles[i] = Triangle(tmpTriangle.A, tmpTriangle.C, tmpTriangle.B)
            i += 1
        # all triangles now have the same orientation, check if it is correct
        tmpVolume = self.Volume()
        if tmpVolume < 0:
            # opposite orientation, flip all triangles
            i = 0
            for tmpTriangle in self.triangles:
                    self.triangles[i] = Triangle(tmpTriangle.A, tmpTriangle.C, tmpTriangle.B)
                    i += 1
        print "OifMesh: Triangulation correct."
        #print strftime("%H:%M:%S", gmtime())
        return 0

    def Surface(self):
        surface = 0.0
        for triangle in self.triangles:
            surface += triangle.Area()
        return surface

    def Volume(self):
        volume = 0.0
        for triangle in self.triangles:
            tmpNormal = GetNTriangle(triangle.A.GetPos(), triangle.B.GetPos(), triangle.C.GetPos())
            tmpNormalLength = Norm(tmpNormal)
            tmpSumZCoords = 1.0 / 3.0 * (triangle.A.GetPos()[2] + triangle.B.GetPos()[2] + triangle.C.GetPos()[2])
            volume -= triangle.Area() * tmpNormal[2] / tmpNormalLength * tmpSumZCoords #tu som to prepisala na minus
        return volume

    def GetNNodes(self):
        return len(self.points)

    def GetNTriangles(self):
        return len(self.triangles)
        
    def GetNEdges(self):
        return len(self.edges)

    def OutputMeshTriangles(self, fileTriangles=None):
        # this is useful after the mesh correction
        # output of mesh nodes can be done from OifCell (this is because their position may change)
        if fileTriangles is None:
            print "OifMesh: No filename provided for triangles. Will use triangles.dat"
            fileTriangles = "triangles.dat"
        outputFile = open(fileTriangles, "w")
        for t in self.triangles:
            outputFile.write(str(t.A.id) + " " + str(t.B.id) + " " + str(t.C.id) + "\n")
        outputFile.close()
        return 0

    def Mirror(self, mirrorX=0, mirrorY=0, mirrorZ=0, outFileName=""):
        if (outFileName == ""):
            print "Cell.Mirror: output meshnodes file for new mesh is missing. Quitting. "
            quit()
        if ((mirrorX!=0 and mirrorX != 1) or (mirrorY!=0 and mirrorY != 1) or (mirrorZ!=0 and mirrorZ != 1)):
            print "Mesh.Mirror: for mirroring only values 0 or 1 are accepted. 1 indicates that the corresponding coordinate will be flipped.  Exiting."
            quit()
        if (mirrorX + mirrorY + mirrorZ > 1):
            print "Mesh.Mirror: flipping allowed only for one axis. Exiting."
            quit()
        if (mirrorX + mirrorY + mirrorZ == 1):
            outFile = open(outFileName, "w")
            for p in self.points:
                coor = p.GetPos()
                if (mirrorX == 1):
                    coor[0] *= -1.0
                if (mirrorY == 1):
                    coor[1] *= -1.0
                if (mirrorZ == 1):
                    coor[2] *= -1.0
                outFile.write(CuStr(coor[0]) + " " + CuStr(coor[1]) + " " + CuStr(coor[2]) + "\n")
            outFile.close()
        return 0


class OifCellType:  # analogous to oif_template
    """ Represents a template for creating elastic objects"""

    def __init__(self, nodesfile="", trianglesfile="", system=None, stretch=(1.0, 1.0, 1.0), ks=0.0, kslin=0.0, kb=0.0, kal=0.0, kag=0.0,
                 kv=0.0, kvisc=0.0, normal=False, checkOrientation=True):
        if system is None:
            print "OifCellType: No system provided. Quitting."
            quit()
        if (nodesfile is "") or (trianglesfile is ""):
            print "OifCellType: One of nodesfile or trianglesfile is missing. Quitting."
            quit()
        if (normal is False):
            print "OifCellType warning: Option normal is not used => membrane collision will not work."
        if (checkOrientation is False):
            print "OifCellType warning: Check of orientation of triangles is switched off."
        if (ks != 0.0) and (kslin != 0.0):
            print "OifCellType: Cannot use linear and nonlinear stretching at the same time. Quitting."
            quit()
        self.system = system
        self.mesh = Mesh(nodesFile=nodesfile, trianglesFile=trianglesfile, system=system, stretch=stretch, normal=normal, checkOrientation=checkOrientation)
        self.localForceInteractions = []
        self.stretch = stretch
        self.ks = ks
        self.kslin = kslin
        self.kb = kb
        self.kal = kal
        self.kag = kag
        self.kv = kv
        self.kvisc = kvisc
        self.normal = normal
        if (ks != 0.0) or (kslin != 0.0) or (kb != 0.0) or (kal != 0.0):
            for angle in self.mesh.angles:
                r0 = Distance(angle.B.GetPos(), angle.C.GetPos())
                phi = AngleBtwTriangles(angle.A.GetPos(), angle.B.GetPos(), angle.C.GetPos(), angle.D.GetPos())
                area1 = AreaTriangle(angle.A.GetPos(), angle.B.GetPos(), angle.C.GetPos())
                area2 = AreaTriangle(angle.D.GetPos(), angle.B.GetPos(), angle.C.GetPos())
                tmpLocalForceInter = OifLocalForces(r0=r0, ks=ks, kslin=kslin, phi0=phi, kb=kb, A01=area1, A02=area2,
                                                    kal=kal, kvisc=kvisc)
                self.localForceInteractions.append([tmpLocalForceInter, [angle.A, angle.B, angle.C, angle.D]])
                self.system.bonded_inter.add(tmpLocalForceInter)
        else:
            print "OifCellType warning: No local interactions created when creating OifCellType."
        if (kag != 0.0) or (kv != 0.0):
            surface = self.mesh.Surface()
            volume = self.mesh.Volume()
            self.globalForceInteraction = OifGlobalForces(A0_g=surface, ka_g=kag, V0=volume, kv=kv)
            self.system.bonded_inter.add(self.globalForceInteraction)
        else:
            print "OifCellType warning: No global interactions created when creating OifCellType."
        self.PrintInfo()

    def PrintInfo(self):
        print "\nThe following OifCellType was created: "
        print "\t nodesFile: " + self.mesh.nodesFile
        print "\t trianglesFile: " + self.mesh.trianglesFile
        print "\t nNodes: " + str(self.mesh.GetNNodes())
        print "\t nTriangles: " + str(self.mesh.GetNTriangles())
        print "\t nEdges: " + str(self.mesh.GetNEdges())
        print "\t ks: " + CuStr(self.ks)
        print "\t kslin: " + CuStr(self.kslin)
        print "\t kb: " + CuStr(self.kb)
        print "\t kal: " + CuStr(self.kal)
        print "\t kag: " + CuStr(self.kag)
        print "\t kvisc: " + CuStr(self.kvisc)
        print "\t normal: " + str(self.normal)
        print "\t stretch: " + str(self.stretch)
        print ""


class OifCell:
    """ Represents a concrete elastic object """

    def __init__(self, cellType=None, origin=None, partType=None, partMass=1.0, rotate=None):
        if cellType is None:
            print "OifCell: No cellType provided. Quitting."
            quit()
        if origin is None:
            print "OifCell: No origin specified. Quitting."
            quit()
        if partType is None:
            print "OifCell: No partType specified. Quitting."
            quit()
        self.cellType = cellType
        self.mesh = cellType.mesh.Copy(origin=origin, partType=partType, partMass=partMass, rotate=rotate)
        self.partMass = partMass
        self.partType = partType
        self.origin = origin
        self.rotate = rotate
        for inter in self.cellType.localForceInteractions:
            espInter = inter[0]
            points = inter[1]
            n_points = len(points)
            if n_points is 2:
                p0 = self.mesh.points[points[0].id]  # Getting PartPoints from id's of FixedPoints
                p1 = self.mesh.points[points[1].id]
                p0.part.add_bond((espInter, p1.partId))
            if n_points is 3:
                p0 = self.mesh.points[points[0].id]
                p1 = self.mesh.points[points[1].id]
                p2 = self.mesh.points[points[2].id]
                p0.part.add_bond((espInter, p1.partId, p2.partId))
            if n_points is 4:
                p0 = self.mesh.points[points[0].id]
                p1 = self.mesh.points[points[1].id]
                p2 = self.mesh.points[points[2].id]
                p3 = self.mesh.points[points[3].id]
                p1.part.add_bond((espInter, p0.partId, p2.partId, p3.partId))
                # Tu treba skontrolovat v akom poradi sa maju zadavat.....

        if ((self.cellType.kag!=0.0) or (self.cellType.kv!=0.0)):
            for triangle in self.mesh.triangles:
                triangle.A.part.add_bond((self.cellType.globalForceInteraction, triangle.B.partId, triangle.C.partId))

        # setting the out_direction interaction for membrane collision
        if self.cellType.mesh.normal is True:
            tmpOutDirectionInteraction = OifOutDirection()
            self.cellType.system.bonded_inter.add(tmpOutDirectionInteraction) # tato interakcia by mohla byt spolocna pre vsetky objekty, ale vytvara sa pre kazdy zvlast
            for p in self.mesh.points:
                p.part.add_bond((tmpOutDirectionInteraction, self.mesh.neighbors[p.id].A.partId, self.mesh.neighbors[p.id].B.partId, self.mesh.neighbors[p.id].C.partId))

        self.PrintInfo()

    def GetOrigin(self):
        center = np.array([0.0, 0.0, 0.0])
        for p in self.mesh.points:
            center += p.GetPos()
        return center/len(self.mesh.points)

    def SetOrigin(self, newOrigin = (0.0, 0.0, 0.0)):
        oldOrigin = self.GetOrigin()
        for p in self.mesh.points:
            newPosition = p.GetPos() - oldOrigin + newOrigin
            p.SetPos(newPosition)

    def GetApproxOrigin(self):
        approxCenter = [0.0, 0.0, 0.0]
        for id in self.mesh.idsExtremalPoints:
            approxCenter += self.mesh.points[id].GetPos()
        return approxCenter/len(self.mesh.idsExtremalPoints)

    def GetOriginFolded(self):
        origin = self.GetOrigin()
        xCoor = np.mod(origin[0], self.cellType.system.box_l[0])
        yCoor = np.mod(origin[1], self.cellType.system.box_l[1])
        zCoor = np.mod(origin[2], self.cellType.system.box_l[2])
        return [xCoor, yCoor, zCoor]

    def GetVelocity(self):
        velocity = np.array([0.0, 0.0, 0.0])
        for p in self.mesh.points:
            velocity += p.GetVel()
        return velocity/len(self.mesh.points)

    def SetVelocity(self, newVelocity = (0.0, 0.0, 0.0)):
        for p in self.mesh.points:
            p.SetVel(newVelocity)

    def PosBounds(self):
        xMin = largeNumber
        xMax = -largeNumber
        yMin = largeNumber
        yMax = -largeNumber
        zMin = largeNumber
        zMax = -largeNumber
        for p in  self.mesh.points:
            coords = p.GetPos()
            if coords[0] < xMin:
               xMin = coords[0]
            if coords[0] > xMax:
                xMax = coords[0]
            if coords[1] < yMin:
                yMin = coords[1]
            if coords[1] > yMax:
                yMax = coords[1]
            if coords[2] < zMin:
                zMin = coords[2]
            if coords[2] > zMax:
                zMax = coords[2]
        return [xMin, xMax, yMin, yMax, zMin, zMax]

    def Surface(self):
        return self.mesh.Surface()

    def Volume(self):
        return self.mesh.Volume()

    def Diameter(self):
        max_distance = 0.0
        nPoints = len(self.mesh.points)
        for i in range(0, nPoints):
            for j in range(i+1, nPoints):
                p1 = self.mesh.points[i].GetPos()
                p2 = self.mesh.points[j].GetPos()
                tmpDist = Distance(p1,p2)
                if (tmpDist > max_distance):
                    max_distance = tmpDist
        return max_distance
        
    def GetNNodes(self):
        return self.mesh.GetNNodes()

    def SetForce(self, newForce = (0.0, 0.0, 0.0)):
        for p in self.mesh.points:
            p.SetForce(newForce)

    def KillMotion(self):
        for p in self.mesh.points:
            p.KillMotion()
            
    def UnkillMotion(self):
        for p in self.mesh.points:
            p.UnkillMotion()

    def OutputVtkPos(self, filename=None):
        if filename is None:
            print "OifCell: No filename provided for vtk output."
            return
        if ".vtk" not in filename:
            print "OifCell warning: A file with vtk format will be written without .vtk extension."
        nPoints = len(self.mesh.points)
        nTriangles = len(self.mesh.triangles)
        outputFile = open(filename, "w")
        outputFile.write("# vtk DataFile Version 3.0\n")
        outputFile.write("Data\n")
        outputFile.write("ASCII\n")
        outputFile.write("DATASET POLYDATA\n")
        outputFile.write("POINTS " + str(nPoints) + " float\n")
        for p in self.mesh.points:
            coords = p.GetPos()
            outputFile.write(CuStr(coords[0]) + " " + CuStr(coords[1]) + " " + CuStr(coords[2]) + "\n")
        outputFile.write("TRIANGLE_STRIPS " + str(nTriangles) + " " + str(4*nTriangles) + "\n")
        for t in self.mesh.triangles:
            outputFile.write("3 " + str(t.A.id) + " " + str(t.B.id) + " " + str(t.C.id) + "\n")
        outputFile.close()

    def OutputVtkPosFolded(self, filename=None):
        if filename is None:
            print "OifCell: No filename provided for vtk output."
            return
        if ".vtk" not in filename:
            print "OifCell warning: A file with vtk format will be written without .vtk extension."
        nPoints = len(self.mesh.points)
        nTriangles = len(self.mesh.triangles)
        
        # get coordinates of the origin
        center = np.array([0.0, 0.0, 0.0])
        for p in self.mesh.points:
            center += p.GetPos()
        center = center/len(self.mesh.points)
        foldX = np.floor(center[0]/self.cellType.system.box_l[0])
        foldY = np.floor(center[1]/self.cellType.system.box_l[1])
        foldZ = np.floor(center[2]/self.cellType.system.box_l[2])
        # these gives how many times is the origin folded in all three directions
        
        outputFile = open(filename, "w")
        outputFile.write("# vtk DataFile Version 3.0\n")
        outputFile.write("Data\n")
        outputFile.write("ASCII\n")
        outputFile.write("DATASET POLYDATA\n")
        outputFile.write("POINTS " + str(nPoints) + " float\n")
        for p in self.mesh.points:
            coords = p.GetPos()
            xCoor = coords[0] - foldX*self.cellType.system.box_l[0]
            yCoor = coords[1] - foldY*self.cellType.system.box_l[1]
            zCoor = coords[2] - foldZ*self.cellType.system.box_l[2]
            outputFile.write(CuStr(xCoor) + " " + CuStr(yCoor) + " " + CuStr(zCoor) + "\n")
        outputFile.write("TRIANGLE_STRIPS " + str(nTriangles) + " " + str(4 * nTriangles) + "\n")
        for t in self.mesh.triangles:
            outputFile.write("3 " + str(t.A.id) + " " + str(t.B.id) + " " + str(t.C.id) + "\n")
        outputFile.close()

    def AppendPointDataToVtk(self, filename=None, dataName=None, data=None, firstAppend=None):
        if filename is None:
            print "OifCell: AppendPointDataToVtk: No filename provided."
            return
        if data is None:
            print "OifCell: AppendPointDataToVtk: No data provided."
            return
        if dataName is None:
            print "OifCell: AppendPointDataToVtk: No dataName provided."
            return
        if firstAppend is None:
            print "OifCell: AppendPointDataToVtk: Need to know whether this is the first data list to be appended for this file."
            return
        nPoints = self.GetNNodes()
        if (len(data) != nPoints):
            print "OifCell: AppendPointDataToVtk: Number of data points does not match number of mesh points."
            return
        outputFile = open(filename, "a")
        if firstAppend is True:
            outputFile.write("POINT_DATA " + str(nPoints) + "\n")
        outputFile.write("SCALARS " + dataName + " float 1\n")
        outputFile.write("LOOKUP_TABLE default\n")
        for p in self.mesh.points:
            outputFile.write(str(data[p.id]) + "\n")
        outputFile.close()

    def OutputRawData(self, filename=None, data=None):
        if filename is None:
            print "OifCell: OutputRawData: No filename provided."
            return
        if data is None:
            print "OifCell: OutputRawData: No data provided."
            return
        nPoints = self.GetNNodes()
        if (len(data) != nPoints):
            print "OifCell: OutputRawData: Number of data points does not match number of mesh points."
            return
        outputFile = open(filename, "w")
        for p in self.mesh.points:
            outputFile.write(" ".join(map(str,data[p.id])) + "\n")
        outputFile.close()

    def OutputMeshPoints(self, filename=None):                   
        if filename is None:
            print "OifCell: No filename provided for mesh nodes output."
            return
        outputFile = open(filename, "w")
        center = self.GetOrigin()
        for p in self.mesh.points:
            coords = p.GetPos() - center
            outputFile.write(CuStr(coords[0]) + " " + CuStr(coords[1]) + " " + CuStr(coords[2]) + "\n")
        outputFile.close()

    def SetMeshPoints(self, filename=None):                   # a toto SetMeshPoints
        if filename is None:
            print "OifCell: No filename provided for SetMeshNodes. "
            return
        center = self.GetOrigin()
        nPoints = self.GetNNodes()

        infile = open(filename, "r")
        nodesCoord = infile.read().split("\n")
        infile.close()
        # removes a blank line at the end of the file if there is any:
        nodesCoord = filter(None, nodesCoord)  # here we have list of lines with triplets of strings
        if len(nodesCoord) is not nPoints:
            print "OifCell: Mesh nodes not set to new positions: number of lines in the file does not equal number of Cell nodes."
            return
        else:
            i = 0
            for line in nodesCoord:  # extracts coordinates from the string line
                line = line.split()
                newPosition = [float(line[0]), float(line[1]), float(line[2])] + center
                self.mesh.points[i].SetPos(newPosition)
                i += 1

    def PrintInfo(self):
        print "\nThe following OifCell was created: "
        print "\t partMass: " + CuStr(self.partMass)
        print "\t partType: " + str(self.partType)
        print "\t rotate: " + str(self.rotate)
        print "\t origin: " + str(self.origin[0]) + " " + str(self.origin[1]) + " " + str(self.origin[2])

    def ElasticForces(self, elasticForces = (0,0,0,0,0,0), fMetric = (0,0,0,0,0,0), vtkFile = None, rawDataFile = None):
        # the order of parameters in elasticForces and in fMetric is as follows (ks, kb, kal, kag, kv, total)
        # vtkFile means that a vtk file for visualisation of elastic forces will be written
        # rawDataFile means that just the elastic forces will be written into the output file

        # chceme tu mat aj kslin???!!!

        for i in range(0,6):
            if (elasticForces[i] != 0) and (elasticForces[i] != 1):
                print "OifCell: ElasticForces: Incorrect argument. elasticForces has to be a sixtuple of 0s and 1s, " \
                      "specifying which elastic forces will be calculated. The order in the sixtuple is (ks, kb, kal, kag, kv, total)."
                return
        for i in range(0,6):
            if (fMetric[i] != 0) and (fMetric[i] != 1):
                print "OifCell: ElasticForces: Incorrect argument. fMetric has to be a sixtuple of 0s and 1s, " \
                      "specifying which fMetric will be calculated. The order in the sixtuple is (ks, kb, kal, kag, kv, total)"
                return
        # calculation of stretching forces and fMetric
        if (elasticForces[0] is 1) or (elasticForces[5] is 1) or (fMetric[0] is 1) or (fMetric[5] is 1):
            # initialize list
            stretchingForcesList = []
            for p in self.mesh.points:
                stretchingForcesList.append([0.0, 0.0, 0.0])
            # calculation uses edges, but results are stored for nodes
            for e in self.mesh.edges:
                aCurrPos = e.A.GetPos()
                bCurrPos = e.B.GetPos()
                aOrigPos = self.cellType.mesh.points[e.A.id].GetPos()
                bOrigPos = self.cellType.mesh.points[e.B.id].GetPos()
                currDist = e.Length()
                origDist = Distance(aOrigPos, bOrigPos)
                tmpStretchingForce = CalcStretchingForce(self.cellType.ks, aCurrPos, bCurrPos, origDist, currDist)
                stretchingForcesList[e.A.id] += tmpStretchingForce
                stretchingForcesList[e.B.id] -= tmpStretchingForce
            # calculation of stretching fMetric, if needed
            if (fMetric[0] is 1):
                ks_fmetric = 0.0
                for p in self.mesh.points:
                    ks_fmetric += Norm(stretchingForcesList[p.id])

        # calculation of bending forces and fMetric
        if (elasticForces[1] is 1) or (elasticForces[5] is 1) or (fMetric[1] is 1) or (fMetric[5] is 1):
            # initialize list
            bendingForcesList = []
            for p in self.mesh.points:
                bendingForcesList.append([0.0, 0.0, 0.0])
            # calculation uses bending incidences, but results are stored for nodes
            for angle in self.mesh.angles:
                aCurrPos = angle.A.GetPos()
                bCurrPos = angle.B.GetPos()
                cCurrPos = angle.C.GetPos()
                dCurrPos = angle.D.GetPos()
                aOrigPos = self.cellType.mesh.points[angle.A.id].GetPos()
                bOrigPos = self.cellType.mesh.points[angle.B.id].GetPos()
                cOrigPos = self.cellType.mesh.points[angle.C.id].GetPos()
                dOrigPos = self.cellType.mesh.points[angle.D.id].GetPos()
                currAngle = angle.Size()
                origAngle = AngleBtwTriangles(aOrigPos, bOrigPos, cOrigPos, dOrigPos)
                tmpBendingForces = CalcBendingForce(self.cellType.kb, aCurrPos, bCurrPos, cCurrPos, dCurrPos, origAngle, currAngle)
                tmpBendingForce1 = np.array([tmpBendingForces[0], tmpBendingForces[1], tmpBendingForces[2]])
                tmpBendingForce2 = np.array([tmpBendingForces[3], tmpBendingForces[4], tmpBendingForces[5]])
                bendingForcesList[angle.A.id] += tmpBendingForce1
                bendingForcesList[angle.B.id] -= 0.5*tmpBendingForce1 + 0.5*tmpBendingForce2
                bendingForcesList[angle.C.id] -= 0.5*tmpBendingForce1 + 0.5*tmpBendingForce2
                bendingForcesList[angle.D.id] += tmpBendingForce2
            # calculation of bending fMetric, if needed
            if (fMetric[1] is 1):
                kb_fmetric = 0.0
                for p in self.mesh.points:
                    kb_fmetric += Norm(bendingForcesList[p.id])

        # calculation of local area forces and fMetric
        if (elasticForces[2] is 1) or (elasticForces[5] is 1) or (fMetric[2] is 1) or (fMetric[5] is 1):
            # initialize list
            localAreaForcesList = []
            for p in self.mesh.points:
                localAreaForcesList.append([0.0, 0.0, 0.0])
            # calculation uses triangles, but results are stored for nodes
            for t in self.mesh.triangles:
                aCurrPos = t.A.GetPos()
                bCurrPos = t.B.GetPos()
                cCurrPos = t.C.GetPos()
                aOrigPos = self.cellType.mesh.points[t.A.id].GetPos()
                bOrigPos = self.cellType.mesh.points[t.B.id].GetPos()
                cOrigPos = self.cellType.mesh.points[t.C.id].GetPos()
                currArea = t.Area()
                origArea = AreaTriangle(aOrigPos, bOrigPos, cOrigPos)
                tmpLocalAreaForces = CalcLocalAreaForce(self.cellType.kal, aCurrPos, bCurrPos, cCurrPos, origArea, currArea)
                localAreaForcesList[t.A.id] += np.array([tmpLocalAreaForces[0], tmpLocalAreaForces[1], tmpLocalAreaForces[2]])
                localAreaForcesList[t.B.id] += np.array([tmpLocalAreaForces[3], tmpLocalAreaForces[4], tmpLocalAreaForces[5]])
                localAreaForcesList[t.C.id] += np.array([tmpLocalAreaForces[6], tmpLocalAreaForces[7], tmpLocalAreaForces[8]])

            # calculation of local area fMetric, if needed
            if (fMetric[2] is 1):
                kal_fmetric = 0.0
                for p in self.mesh.points:
                    kal_fmetric += Norm(localAreaForcesList[p.id])

        # calculation of global area forces and fMetric
        if (elasticForces[3] is 1) or (elasticForces[5] is 1) or (fMetric[3] is 1) or (fMetric[5] is 1):
            # initialize list
            globalAreaForcesList = []
            for p in self.mesh.points:
                globalAreaForcesList.append([0.0, 0.0, 0.0])
            # calculation uses triangles, but results are stored for nodes
            for t in self.mesh.triangles:
                aCurrPos = t.A.GetPos()
                bCurrPos = t.B.GetPos()
                cCurrPos = t.C.GetPos()
                currSurface = self.mesh.Surface()
                origSurface = self.cellType.mesh.Surface()
                tmpGlobalAreaForces = CalcGlobalAreaForce(self.cellType.kag, aCurrPos, bCurrPos, cCurrPos, origSurface, currSurface)
                globalAreaForcesList[t.A.id] += np.array([tmpGlobalAreaForces[0], tmpGlobalAreaForces[1], tmpGlobalAreaForces[2]])
                globalAreaForcesList[t.B.id] += np.array([tmpGlobalAreaForces[3], tmpGlobalAreaForces[4], tmpGlobalAreaForces[5]])
                globalAreaForcesList[t.C.id] += np.array([tmpGlobalAreaForces[6], tmpGlobalAreaForces[7], tmpGlobalAreaForces[8]])
            # calculation of global area fMetric, if needed
            if (fMetric[3] is 1):
                kag_fmetric = 0.0
                for p in self.mesh.points:
                    kag_fmetric += Norm(globalAreaForcesList[p.id])

        # calculation of volume forces and fMetric
        if (elasticForces[4] is 1) or (elasticForces[5] is 1) or (fMetric[4] is 1) or (fMetric[5] is 1):
            # initialize list
            volumeForcesList = []
            for p in self.mesh.points:
                volumeForcesList.append([0.0, 0.0, 0.0])
            # calculation uses triangles, but results are stored for nodes
            for t in self.mesh.triangles:
                aCurrPos = t.A.GetPos()
                bCurrPos = t.B.GetPos()
                cCurrPos = t.C.GetPos()
                currVolume = self.mesh.Volume()
                origVolume = self.cellType.mesh.Volume()
                tmpVolumeForce = CalcVolumeForce(self.cellType.kv, aCurrPos, bCurrPos, cCurrPos, origVolume, currVolume)
                volumeForcesList[t.A.id] += tmpVolumeForce
                volumeForcesList[t.B.id] += tmpVolumeForce
                volumeForcesList[t.C.id] += tmpVolumeForce
            # calculation of volume fMetric, if needed
            if (fMetric[4] is 1):
                kv_fmetric = 0.0
                for p in self.mesh.points:
                    kv_fmetric += Norm(volumeForcesList[p.id])

        # calculation of total elastic forces and fMetric
        if (elasticForces[5] is 1) or (fMetric[5] is 1):
            elasticForcesList = []
            for p in self.mesh.points:
                totalElasticForces = stretchingForcesList[p.id] + bendingForcesList[p.id] + localAreaForcesList[p.id] \
                                     + globalAreaForcesList[p.id] + volumeForcesList[p.id]
                elasticForcesList.append(totalElasticForces)
            # calculation of total fMetric, if needed
            if (fMetric[5] is 1):
                total_fmetric = 0.0
                for p in self.mesh.points:
                    total_fmetric += Norm(elasticForcesList[p.id])

        # calculate norms of resulting forces
        if (elasticForces[0] + elasticForces[1] + elasticForces[2] + elasticForces[3] + elasticForces[4] + elasticForces[5]) is not 0:
            if (elasticForces[0] is 1):
                stretchingForcesNormsList = []
                for p in self.mesh.points:
                    stretchingForcesNormsList.append(Norm(stretchingForcesList[p.id]))
            if (elasticForces[1] is 1):
                bendingForcesNormsList = []
                for p in self.mesh.points:
                    bendingForcesNormsList.append(Norm(bendingForcesList[p.id]))
            if (elasticForces[2] is 1):
                localAreaForcesNormsList = []
                for p in self.mesh.points:
                    localAreaForcesNormsList.append(Norm(localAreaForcesList[p.id]))
            if (elasticForces[3] is 1):
                globalAreaForcesNormsList = []
                for p in self.mesh.points:
                    globalAreaForcesNormsList.append(Norm(globalAreaForcesList[p.id]))
            if (elasticForces[4] is 1):
                volumeForcesNormsList = []
                for p in self.mesh.points:
                    volumeForcesNormsList.append(Norm(volumeForcesList[p.id]))
            if (elasticForces[5] is 1):
                elasticForcesNormsList = []
                for p in self.mesh.points:
                    elasticForcesNormsList.append(Norm(elasticForcesList[p.id]))

        # output vtk (folded)
        if vtkFile is not None:
            if (elasticForces == (0,0,0,0,0,0)):
                print "OifCell: ElasticForces: The option elasticForces was not used. " \
                      "Nothing to output to vtk file."
                return
            self.OutputVtkPosFolded(vtkFile)
            first = True
            if (elasticForces[0] is 1):
                self.AppendPointDataToVtk(filename=vtkFile, dataName="ks_fMetric", data=stretchingForcesNormsList, firstAppend=first)
                first = False
            if (elasticForces[1] is 1):
                self.AppendPointDataToVtk(filename=vtkFile, dataName="kb_fMetric", data=bendingForcesNormsList, firstAppend=first)
                first = False
            if (elasticForces[2] is 1):
                self.AppendPointDataToVtk(filename=vtkFile, dataName="kal_fMetric", data=localAreaForcesNormsList, firstAppend=first)
                first = False
            if (elasticForces[3] is 1):
                self.AppendPointDataToVtk(filename=vtkFile, dataName="kag_fMetric", data=globalAreaForcesNormsList, firstAppend=first)
                first = False
            if (elasticForces[4] is 1):
                self.AppendPointDataToVtk(filename=vtkFile, dataName="kav_fMetric", data=volumeForcesNormsList, firstAppend=first)
                first = False
            if (elasticForces[5] is 1):
                self.AppendPointDataToVtk(filename=vtkFile, dataName="total_fMetric", data=elasticForcesNormsList, firstAppend=first)
                first = False

        # output raw data
        if rawDataFile is not None:
            if (elasticForces[0] + elasticForces[1] + elasticForces[2] + elasticForces[3] + elasticForces[4] + elasticForces[5]) is not 1:
                print "OifCell: ElasticForces: Only one type of elastic forces can be written into one rawDataFile. " \
                      "If you need several, please call OifCell.ElasticForces multiple times - once per elastic force."
                return
            if (elasticForces[0] is 1):
                self.OutputRawData(filename=rawDataFile, data=stretchingForcesList)
            if (elasticForces[1] is 1):
                self.OutputRawData(filename=rawDataFile, data=bendingForcesList)
            if (elasticForces[2] is 1):
                self.OutputRawData(filename=rawDataFile, data=localAreaForcesList)
            if (elasticForces[3] is 1):
                self.OutputRawData(filename=rawDataFile, data=globalAreaForcesList)
            if (elasticForces[4] is 1):
                self.OutputRawData(filename=rawDataFile, data=volumeForcesList)
            if (elasticForces[5] is 1):
                self.OutputRawData(filename=rawDataFile, data=elasticForcesList)

        # return fMetric
        if fMetric[0] + fMetric[1] + fMetric[2] + fMetric[3] + fMetric[4] + fMetric[5] > 0:
            results = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            if fMetric[0] is 1:
                results[0] = ks_fmetric
            if fMetric[1] is 1:
                results[1] = kb_fmetric
            if fMetric[2] is 1:
                results[2] = kal_fmetric
            if fMetric[3] is 1:
                results[3] = kag_fmetric
            if fMetric[4] is 1:
                results[4] = kv_fmetric
            if fMetric[5] is 1:
                results[5] = total_fmetric
            return results
        else:
            return 0
