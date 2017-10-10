import numpy as np
import math

smallEpsilon = 0.000000001
largeNumber = 10000000.0
pi = 3.1415926535897931
outputPrecision = 14


def CuStr(realn):
	return str('{:.{prec}f}'.format(realn, prec=outputPrecision))


def GetNTriangle(a, b, c):
    n = [0.0, 0.0, 0.0]
    n[0] = (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
    n[1] = (b[2] - a[2]) * (c[0] - a[0]) - (b[0] - a[0]) * (c[2] - a[2])
    n[2] = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    return n

def Norm(v):
    norm = np.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
    return norm

def Distance(a, b):
    dist = np.sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) + (a[2] - b[2]) * (a[2] - b[2]))
    return dist

def AreaTriangle(a, b, c):
    n = GetNTriangle(a, b, c)
    area = 0.5 * Norm(n)
    return area

def AngleBtwTriangles(P1, P2, P3, P4):
    n1 = GetNTriangle(P2, P1, P3)
    n2 = GetNTriangle(P2, P3, P4)
    tmp11 = n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2]
    tmp11 = tmp11 * abs(tmp11)

    tmp22 = n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2]
    tmp33 = n2[0] * n2[0] + n2[1] * n2[1] + n2[2] * n2[2]
    tmp11 = tmp11 / (tmp22 * tmp33)

    if tmp11 > 0:
        tmp11 = np.sqrt(tmp11)
    else:
        tmp11 = - np.sqrt(- tmp11)

    if tmp11 >= 1.0:
        tmp11 = 0.0
    elif tmp11 <= -1.:
        tmp11 = pi

    phi = pi - math.acos(tmp11)

    tmp11 = -(n1[0] * P1[0] + n1[1] * P1[1] + n1[2] * P1[2])
    if n1[0] * P4[0] + n1[1] * P4[1] + n1[2] * P4[2] + tmp11 < 0:
        phi = 2.0 * pi - phi
    return phi

def DiscardEpsilon(x):
    if (x > -smallEpsilon and x < smallEpsilon):
        res = 0.0
    else:
        res = x
    return res

def KS(lambd):
    # Defined by (19) from Dupin2007
    res = (pow(lambd, 0.5) + pow(lambd, -2.5)) / (lambd + pow(lambd, -3.))
    return res

def CalcStretchingForce(ks, pA, pB, dist0, dist):
    # this has to correspond to the calculation in oif_local_forces.hpp: calc_oif_local
    # as of now, corresponds to git commit f156f9b44dcfd3cef9dd5537a1adfc903ac4772a
    dr = dist - dist0
    # nonlinear stretching:
    lambd = 1.0 * dist / dist0
    fac = ks * KS(lambd) * dr       # no negative sign here! different from C implementation due to reverse order of vector subtraction
    f = fac * (np.array(pB) - np.array(pA)) / dist
    return f

def CalcLinearStretchingForce(ks, pA, pB, dist0, dist):
    dr = dist - dist0
    fac = ks * dr                   # no negative sign here! different from C implementation due to reverse order of vector subtraction
    f = fac * (np.array(pB) - np.array(pA)) / dist
    return f

def CalcBendingForce(kb, pA, pB, pC, pD, phi0, phi):
    # this has to correspond to the calculation in oif_local_forces.hpp: calc_oif_local
    # as of now, corresponds to git commit f156f9b44dcfd3cef9dd5537a1adfc903ac4772a
    n1 = GetNTriangle(pB, pA, pC)
    dn1 = Norm(n1)

    n2 = GetNTriangle(pB, pC, pD)
    dn2 = Norm(n2)

    angles = (phi - phi0) / phi0
    fac = kb * angles

    f1 = fac * np.array(n1) / dn1
    f2 = fac * np.array(n2) / dn2
    f = [f1[0], f1[1], f1[2], f2[0], f2[1], f2[2]]
    return f

def CalcLocalAreaForce(kal, pA, pB, pC, A0, A):
    # this has to correspond to the calculation in oif_local_forces.hpp: calc_oif_local
    # except for division by 3 - each triangle enters this calculation once, while each triangle enters the
    # calc_oif_local three times
    # as of now, corresponds to git commit f156f9b44dcfd3cef9dd5537a1adfc903ac4772a

    centroid = np.array((pA + pB + pC) / 3.0)
    deltaA = A - A0
    ta = centroid - pA
    taNorm = Norm(ta)
    tb = centroid - pB
    tbNorm = Norm(tb)
    tc = centroid - pC
    tcNorm = Norm(tc)
    commonFactor = kal * deltaA / (taNorm * taNorm + tbNorm * tbNorm + tcNorm * tcNorm)

    # local area force for first node
    f1 = commonFactor * ta

    # local area force for second node
    f2 = commonFactor * tb

    # local area force for third node
    f3 = commonFactor * tc

    f = [f1[0], f1[1], f1[2], f2[0], f2[1], f2[2], f3[0], f3[1], f3[2]]
    return f

def CalcGlobalAreaForce(kag, pA, pB, pC, Ag0, Ag):
    # this has to correspond to the calculation in oif_global_forces.hpp: add_oif_global_forces
    # as of now, corresponds to git commit f156f9b44dcfd3cef9dd5537a1adfc903ac4772a
    centroid = np.array((pA + pB + pC) / 3.0)
    delta = Ag - Ag0
    ta = centroid - pA
    taNorm = Norm(ta)
    tb = centroid - pB
    tbNorm = Norm(tb)
    tc = centroid - pC
    tcNorm = Norm(tc)
    A = AreaTriangle(pA, pB, pC)
    commonFactor = kag * A * delta / (taNorm * taNorm + tbNorm * tbNorm + tcNorm * tcNorm)
    #print "fac from py: "+str(commonFactor)

    # global area force for first node
    f1 = commonFactor * ta

    # global area force for second node
    f2 = commonFactor * tb

    # global area force for third node
    f3 = commonFactor * tc

    f = [f1[0], f1[1], f1[2], f2[0], f2[1], f2[2], f3[0], f3[1], f3[2]]
    return f

def CalcVolumeForce(kv, pA, pB, pC, V0, V):
    # this has to correspond to the calculation in oif_global_forces.hpp: add_oif_global_forces
    # as of now, corresponds to git commit f156f9b44dcfd3cef9dd5537a1adfc903ac4772a
    n = GetNTriangle(pA, pB, pC)
    dn = Norm(n)
    vv = (V - V0) / V0
    A = AreaTriangle(pA, pB, pC)
    f = kv * vv * A * np.array(n) / (dn * 3.0)
    return f

def OutputVtkRhomboid(corner, a, b, c, outFile):
    if ".vtk" not in outFile:
        print "OutputVtkRhomboid warning: A file with vtk format will be written without .vtk extension."
    outputFile = open(outFile, "w")
    outputFile.write("# vtk DataFile Version 3.0\n")
    outputFile.write("Data\n")
    outputFile.write("ASCII\n")
    outputFile.write("DATASET POLYDATA\n")
    outputFile.write("POINTS 8 float\n")

    outputFile.write(str(corner[0]) + " " + str(corner[1]) + " " + str(corner[2]) + "\n")
    outputFile.write(str(corner[0] + a[0]) + " " + str(corner[1] + a[1]) + " " + str(corner[2] + a[2]) + "\n")
    outputFile.write(str(corner[0] + a[0] + b[0]) + " " + str(corner[1] + a[1] + b[1]) + " " + str(corner[2] + a[2] + b[2]) + "\n")
    outputFile.write(str(corner[0] + b[0]) + " " + str(corner[1] + b[1]) + " " + str(corner[2] + b[2]) + "\n")

    outputFile.write(str(corner[0] + c[0]) + " " + str(corner[1] + c[1]) + " " + str(corner[2] + c[2]) + "\n")
    outputFile.write(str(corner[0] + a[0] + c[0]) + " " + str(corner[1] + a[1] + c[1]) + " " + str(corner[2] + a[2] + c[2]) + "\n")
    outputFile.write(str(corner[0] + a[0] + b[0] + c[0]) + " " + str(corner[1] + a[1] + b[1] + c[1]) + " " + str(corner[2] + a[2] + b[2] + c[2]) + "\n")
    outputFile.write(str(corner[0] + b[0] + c[0]) + " " + str(corner[1] + b[1] + c[1]) + " " + str(corner[2] + b[2] + c[2]) + "\n")

    outputFile.write("POLYGONS 6 30\n")
    outputFile.write("4 0 1 2 3\n")
    outputFile.write("4 4 5 6 7\n")
    outputFile.write("4 0 1 5 4\n")
    outputFile.write("4 2 3 7 6\n")
    outputFile.write("4 0 4 7 3\n")
    outputFile.write("4 1 2 6 5")
    outputFile.close()
    return 0

def OutputVtkCylinder(center, axis, length, radius, n, outFile):
    # L is the half height of the cylinder
    # only vertical cylinders are supported for now, i.e. with normal (0.0, 0.0, 1.0)

    if ".vtk" not in outFile:
        print "OutputVtkCylinder warning: A file with vtk format will be written without .vtk extension."
    checkAxis = True
    if axis[0]!=0.0:
        checkAxis = False
    if axis[1]!=0.0:
        checkAxis = False
    if axis[2]==0.0:
        checkAxis = False
    if checkAxis is False:
        print "OutputVtkCylinder: Output for this type of cylinder is not supported yet."
        return
    axisZ = 1.0

    # setting points on perimeter
    alpha = 2 * pi / n
    points = 2 * n

    # shift center to the bottom circle
    p1 = center - length * np.array(axis)

    outputFile = open(outFile, "w")
    outputFile.write("# vtk DataFile Version 3.0\n")
    outputFile.write("Data\n")
    outputFile.write("ASCII\n")
    outputFile.write("DATASET POLYDATA\n")
    outputFile.write("POINTS " + str(points) + " float\n")
    for i in range(0,n):
        outputFile.write(str(p1[0] + radius*np.cos(i*alpha)) + " " + str(p1[1] + radius*np.sin(i*alpha)) + " " + str(p1[2]) + "\n")
    for i in range(0,n):
        outputFile.write(str(p1[0] + radius*np.cos(i*alpha)) + " " + str(p1[1] + radius*np.sin(i*alpha)) + " " + str(p1[2]+2*length*axisZ) + "\n")
    outputFile.write("POLYGONS " + str(n+2) + " " + str(5*n+(n+1)*2) + "\n")

    # writing bottom "circle"
    outputFile.write(str(n) + " ")
    for i in range(0, n - 1):
        outputFile.write(str(i) + " ")
    outputFile.write(str(n-1) + "\n")

    # writing top "circle"
    outputFile.write(str(n) + " ")
    for i in range(0, n - 1):
        outputFile.write(str(i+n) + " ")
    outputFile.write(str(2*n - 1) + "\n")

    # writing sides - rectangles
    for i in range(0, n - 1):
        outputFile.write("4 " + str(i) + " " + str(i+1) + " " + str(i+n+1) + " " + str(i+n) + "\n")
    outputFile.write("4 " + str(n-1) + " " + str(0) + " " + str(n) + " " + str(2*n-1) + "\n")

    outputFile.close()
    return 0

def OutputVtkLines(lines, outFile):
    # lines is a list of pairs of points p1, p2
    # each pair represents a line segment to output to vtk
    # each line in lines contains 6 floats: p1x, p1y, p1z, p2x, p2y, p2z

    if ".vtk" not in outFile:
        print "OutputVtkLines warning: A file with vtk format will be written without .vtk extension."

    nlines = len(lines)

    outputFile = open(outFile, "w")
    outputFile.write("# vtk DataFile Version 3.0\n")
    outputFile.write("Data\n")
    outputFile.write("ASCII\n")
    outputFile.write("DATASET POLYDATA\n")
    outputFile.write("POINTS " + str(2*nlines) + " float\n")
    for i in range(0,nlines):
        oneline = lines[i]
        outputFile.write(str(oneline[0]) + " " + str(oneline[1]) + " " + str(oneline[2]) + "\n")
        outputFile.write(str(oneline[3]) + " " + str(oneline[4]) + " " + str(oneline[5]) + "\n")
    outputFile.write("LINES " + str(nlines) + " " + str(3 *nlines) + "\n")
    for i in range(0,nlines):
        outputFile.write(str(2) + " " + str(2*i) + " " + str(2*i+1) + "\n")

    outputFile.close()
    return 0

def GetInterpolatedVelocity(position,lbf):
    # position is a vector [x,y,z] in which we want the fluid velocity
    # WARNING: this is only for agrid = 1

    xLower = math.floor(position[0])
    xUpper = math.ceil(position[0])
    yLower = math.floor(position[1])
    yUpper = math.ceil(position[1])
    zLower = math.floor(position[2])
    zUpper = math.ceil(position[2])

    vel000 = np.array(lbf[xLower, yLower, zLower].velocity)
    vel100 = np.array(lbf[xUpper, yLower, zLower].velocity)
    vel010 = np.array(lbf[xLower, yUpper, zLower].velocity)
    vel001 = np.array(lbf[xLower, yLower, zUpper].velocity)
    vel110 = np.array(lbf[xUpper, yUpper, zLower].velocity)
    vel101 = np.array(lbf[xUpper, yLower, zUpper].velocity)
    vel011 = np.array(lbf[xLower, yUpper, zUpper].velocity)
    vel111 = np.array(lbf[xUpper, yUpper, zUpper].velocity)

    x = position[0] - xLower
    y = position[1] - yLower
    z = position[2] - zLower

    interpVel = vel000 * (1 - x) * (1 - y) * (1 - z) \
        + vel001 * (1 - x) * (1 - y) * z \
        + vel010 * (1 - x) * y * (1 - z) \
        + vel100 * x * (1 - y) * (1 - z) \
        + vel110 * x * y * (1 - z) \
        + vel101 * x * (1 - y) * z \
        + vel011 * (1 - x) * y * z \
        + vel111 * x * y * z

    return interpVel
