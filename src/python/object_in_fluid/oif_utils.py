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
import numpy as np
import math
from espressomd.shapes import Rhomboid
from espressomd.shapes import Cylinder

small_epsilon = 0.000000001
large_number = 10000000.0
output_precision = 14


def custom_str(realn):
    return str('{:.{prec}f}'.format(realn, prec=output_precision))


def get_triangle_normal(a, b, c):
    """
    Returns the normal vector of a triangle given by points a,b,c.

    Parameters
    ----------
    a : (3,) array_like of :obj:`float`
          Point a
    b : (3,) array_like of :obj:`float`
          Point b
    c : (3,) array_like of :obj:`float`
          Point c
    """
    n = [0.0, 0.0, 0.0]
    n[0] = (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
    n[1] = (b[2] - a[2]) * (c[0] - a[0]) - (b[0] - a[0]) * (c[2] - a[2])
    n[2] = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    return np.array(n)


def norm(vect):
    """
    Returns the norm of a vector.

    Parameters
    ----------
    vect : (3,) array_like of :obj:`float`
          Input vector

    """
    v = np.array(vect)
    return np.sqrt(np.dot(v, v))


def vec_distance(a, b):
    """
    Returns the length of vector between points a and b.

    Parameters
    ----------
    a : (3,) array_like of :obj:`float`
          Point a
    b : (3,) array_like of :obj:`float`
          Point b

    """
    return norm(np.array(a) - np.array(b))


def area_triangle(a, b, c):
    """
    Returns the area of a triangle given by points a, b, c.


    Parameters
    ----------
    a : (3,) array_like of :obj:`float`
          Point a
    b : (3,) array_like of :obj:`float`
          Point b
    c : (3,) array_like of :obj:`float`
          Point c

    """
    n = get_triangle_normal(a, b, c)
    area = 0.5 * norm(n)
    return area


def angle_btw_triangles(P1, P2, P3, P4):
    """
    Returns the size of an angle between triangles given by points P2, P1, P3 and P2, P3, P4.

    Parameters
    ----------
    P1 : (3,) array_like of :obj:`float`
          Point P1
    P2 : (3,) array_like of :obj:`float`
          Point P2
    P3 : (3,) array_like of :obj:`float`
          Point P3
    P4 : (3,) array_like of :obj:`float`
          Point P4

    """
    n1 = get_triangle_normal(P2, P1, P3)
    n2 = get_triangle_normal(P2, P3, P4)
    tmp11 = np.dot(n1, n2) / (np.linalg.norm(n1) * np.linalg.norm(n2))

    if tmp11 >= 1.0:
        tmp11 = 1.0
    elif tmp11 <= -1.0:
        tmp11 = -1.0

    phi = np.pi - math.acos(tmp11)

    if (np.dot(n1, np.array(P4)) - np.dot(n1, np.array(P1))) < 0:
        phi = 2.0 * np.pi - phi
    return phi


def discard_epsilon(x):
    """
    Returns zero if the argument is too small.

    Parameters
    ----------
    x : :obj:`float`
          real number

    """
    if (x > -small_epsilon and x < small_epsilon):
        res = 0.0
    else:
        res = x
    return res


def oif_neo_hookean_nonlin(lambd):
    """
    Defines NeoHookean nonlinearity.

    Parameters
    ----------
    lambd : :obj:`float`
          real number
    """
    # Defined by (19) from Dupin2007
    res = (pow(lambd, 0.5) + pow(lambd, -2.5)) / (lambd + pow(lambd, -3.))
    return res


def oif_calc_stretching_force(ks, pA, pB, dist0, dist):
    """
    Calculates nonlinear stretching forces between two points on an edge.

    Parameters
    ----------
    ks : :obj:`float`
          coefficient of the stretching, spring stiffness
    pA : (3,) array_like of :obj:`float`
          position of the first particle
    pB : (3,) array_like of :obj:`float`
          position of the second particle
    dist0 : :obj:`float`
          relaxed distance btw particles
    dist : :obj:`float`
          current distance btw particles


    """
    # this has to correspond to the calculation in oif_local_forces.hpp: calc_oif_local
    # as of now, corresponds to git commit
    # f156f9b44dcfd3cef9dd5537a1adfc903ac4772a
    dr = dist - dist0
    # nonlinear stretching:
    lambd = 1.0 * dist / dist0
    fac = ks * oif_neo_hookean_nonlin(lambd) * dr
    # no negative sign here! different from C implementation
    # due to reverse order of vector subtraction
    f = fac * (np.array(pB) - np.array(pA)) / dist
    return f


def oif_calc_linear_stretching_force(ks, pA, pB, dist0, dist):
    """
    Calculates linear stretching forces between two points on an edge.

    Parameters
    ----------
    ks : :obj:`float`
          coefficient of the stretching, spring stiffness
    pA : (3,) array_like of :obj:`float`
          position of the first particle
    pB : (3,) array_like of :obj:`float`
          position of the second particle
    dist0 : :obj:`float`
          relaxed distance btw particles
    dist : :obj:`float`
          current distance btw particles

    """
    dr = dist - dist0
    fac = ks * dr
    # no negative sign here! different from C implementation due to
    # reverse order of vector subtraction
    f = fac * (np.array(pB) - np.array(pA)) / dist
    return f


def oif_calc_bending_force(kb, pA, pB, pC, pD, phi0, phi):
    """
    Calculates bending forces for four points on two adjacent triangles.

    Parameters
    ----------
    kb : :obj:`float`
          coefficient of the stretching, spring stiffness
    pA : (3,) array_like of :obj:`float`
          position of the first particle
    pB : (3,) array_like of :obj:`float`
          position of the second particle
    pC : (3,) array_like of :obj:`float`
          position of the third particle
    pD : (3,) array_like of :obj:`float`
          position of the fourth particle
    phi0 : :obj:`float`
          relaxed angle btw two triangles
    phi : :obj:`float`
          current angle btw two triangles

    """
    # this has to correspond to the calculation in oif_local_forces.hpp: calc_oif_local
    # as of now, corresponds to git commit
    # f156f9b44dcfd3cef9dd5537a1adfc903ac4772a
    n1 = get_triangle_normal(pB, pA, pC)
    n2 = get_triangle_normal(pB, pC, pD)

    angles = (phi - phi0) / phi0
    fac = kb * angles

    f1 = fac * np.array(n1) / norm(n1)
    f2 = fac * np.array(n2) / norm(n2)
    f = [f1[0], f1[1], f1[2], f2[0], f2[1], f2[2]]
    return f


def oif_calc_local_area_force(kal, pA, pB, pC, A0, A):
    """
    Calculates local area forces between three points in one triangle.

    Parameters
    ----------
    kal : :obj:`float`
          coefficient of the stretching, spring stiffness
    pA : (3,) array_like of :obj:`float`
          position of the first particle
    pB : (3,) array_like of :obj:`float`
          position of the second particle
    pC : (3,) array_like of :obj:`float`
          position of the third particle
    A0 : :obj:`float`
          relaxed area of the triangle
    A : :obj:`float`
          current area of the triangle

    """
    # this has to correspond to the calculation in oif_local_forces.hpp: calc_oif_local
    # except for division by 3 - each triangle enters this calculation once, while each triangle enters the
    # calc_oif_local three times
    # as of now, corresponds to git commit
    # f156f9b44dcfd3cef9dd5537a1adfc903ac4772a

    centroid = np.array((pA + pB + pC) / 3.0)
    delta_area = A - A0
    ta = centroid - pA
    ta_norm = norm(ta)
    tb = centroid - pB
    tb_norm = norm(tb)
    tc = centroid - pC
    tc_norm = norm(tc)
    common_factor = kal * delta_area / \
        (ta_norm * ta_norm + tb_norm * tb_norm + tc_norm * tc_norm)

    # local area force for first node
    f1 = common_factor * ta

    # local area force for second node
    f2 = common_factor * tb

    # local area force for third node
    f3 = common_factor * tc

    f = [f1[0], f1[1], f1[2], f2[0], f2[1], f2[2], f3[0], f3[1], f3[2]]
    return f


def oif_calc_global_area_force(kag, pA, pB, pC, Ag0, Ag):
    """
    Calculates global area forces between three points in a triangle.

    Parameters
    ----------
    kag : :obj:`float`
          coefficient of the stretching, spring stiffness
    pA : (3,) array_like of :obj:`float`
          position of the first particle
    pB : (3,) array_like of :obj:`float`
          position of the second particle
    pC : (3,) array_like of :obj:`float`
          position of the third particle
    Ag0 : :obj:`float`
          relaxed surface area of the cell
    Ag : :obj:`float`
          current surface area of the cell

    """
    # this has to correspond to the calculation in oif_global_forces.hpp: add_oif_global_forces
    # as of now, corresponds to git commit
    # f156f9b44dcfd3cef9dd5537a1adfc903ac4772a
    centroid = np.array((pA + pB + pC) / 3.0)
    delta = Ag - Ag0
    ta = centroid - pA
    ta_norm = norm(ta)
    tb = centroid - pB
    tb_norm = norm(tb)
    tc = centroid - pC
    tc_norm = norm(tc)
    A = area_triangle(pA, pB, pC)
    common_factor = kag * A * delta / \
        (ta_norm * ta_norm + tb_norm * tb_norm + tc_norm * tc_norm)

    # global area force for first node
    f1 = common_factor * ta

    # global area force for second node
    f2 = common_factor * tb

    # global area force for third node
    f3 = common_factor * tc

    f = [f1[0], f1[1], f1[2], f2[0], f2[1], f2[2], f3[0], f3[1], f3[2]]
    return f


def oif_calc_volume_force(kv, pA, pB, pC, V0, V):
    """
    Calculates volume forces for three points in a triangle.

    Parameters
    ----------
    kv : :obj:`float`
          coefficient of the stretching, spring stiffness
    pA : (3,) array_like of :obj:`float`
          position of the first particle
    pB : (3,) array_like of :obj:`float`
          position of the second particle
    pC : (3,) array_like of :obj:`float`
          position of the third particle
    V0 : :obj:`float`
          relaxed volume of the cell
    V : :obj:`float`
          current volume of the cell

    """
    # this has to correspond to the calculation in oif_global_forces.hpp: add_oif_global_forces
    # as of now, corresponds to git commit
    # f156f9b44dcfd3cef9dd5537a1adfc903ac4772a
    n = get_triangle_normal(pA, pB, pC)
    dn = norm(n)
    vv = (V - V0) / V0
    A = area_triangle(pA, pB, pC)
    f = kv * vv * A * np.array(n) / (dn * 3.0)
    return f


def output_vtk_rhomboid(rhom_shape, out_file):
    """
    Outputs the VTK files for visualisation of a rhomboid in e.g. Paraview.

    Parameters
    ----------
    rhom_shape : :obj:`shape`
          rhomboid shape
    out_file : :obj:`str`
          filename for the output

    """

    corner = rhom_shape.corner
    a = rhom_shape.a
    b = rhom_shape.b
    c = rhom_shape.c

    output_file = open(out_file, "w")
    output_file.write("# vtk DataFile Version 3.0\n")
    output_file.write("Data\n")
    output_file.write("ASCII\n")
    output_file.write("DATASET POLYDATA\n")
    output_file.write("POINTS 8 float\n")

    output_file.write(str(corner[0]) + " " + str(
        corner[1]) + " " + str(corner[2]) + "\n")
    output_file.write(str(corner[0] + a[0]) + " " + str(
        corner[1] + a[1]) + " " + str(corner[2] + a[2]) + "\n")
    output_file.write(str(corner[0] + a[0] + b[0]) + " " + str(corner[1] + a[1] + b[1]) + " " +
                      str(corner[2] + a[2] + b[2]) + "\n")
    output_file.write(str(corner[0] + b[0]) + " " + str(
        corner[1] + b[1]) + " " + str(corner[2] + b[2]) + "\n")

    output_file.write(str(corner[0] + c[0]) + " " + str(
        corner[1] + c[1]) + " " + str(corner[2] + c[2]) + "\n")
    output_file.write(str(corner[0] + a[0] + c[0]) + " " + str(corner[1] + a[1] + c[1]) + " " +
                      str(corner[2] + a[2] + c[2]) + "\n")
    output_file.write(str(corner[0] + a[0] + b[0] + c[0]) + " " + str(corner[1] + a[1] + b[1] + c[1]) + " " +
                      str(corner[2] + a[2] + b[2] + c[2]) + "\n")
    output_file.write(str(corner[0] + b[0] + c[0]) + " " + str(corner[1] + b[1] + c[1]) + " " +
                      str(corner[2] + b[2] + c[2]) + "\n")

    output_file.write("POLYGONS 6 30\n")
    output_file.write("4 0 1 2 3\n")
    output_file.write("4 4 5 6 7\n")
    output_file.write("4 0 1 5 4\n")
    output_file.write("4 2 3 7 6\n")
    output_file.write("4 0 4 7 3\n")
    output_file.write("4 1 2 6 5")
    output_file.close()
    return 0


def output_vtk_cylinder(cyl_shape, n, out_file):
    """
    Outputs the VTK files for visualisation of a cylinder in e.g. Paraview.

    Parameters
    ----------
    cyl_shape : :obj:`shape`
          cylindrical shape
    n : :obj:`int`
          number of discretization sections
    out_file : :obj:`str`
          filename for the output

    """
    # length is the full height of the cylinder (note: used to be just half in the previous versions)
    # only vertical cylinders are supported for now, i.e. with normal (0.0,
    # 0.0, 1.0)

    axis = cyl_shape.axis
    length = cyl_shape.length
    radius = cyl_shape.radius
    center = cyl_shape.center

    check_axis = True
    if axis[0] != 0.0:
        check_axis = False
    if axis[1] != 0.0:
        check_axis = False
    if axis[2] == 0.0:
        check_axis = False
    if check_axis is False:
        raise Exception(
            "output_vtk_cylinder: Output for this type of cylinder is not supported yet.")
    axisZ = 1.0

    # setting points on perimeter
    alpha = 2 * np.pi / n
    points = 2 * n

    # shift center to the bottom circle
    p1 = center - length * np.array(axis) / 2.0

    output_file = open(out_file, "w")
    output_file.write("# vtk DataFile Version 3.0\n")
    output_file.write("Data\n")
    output_file.write("ASCII\n")
    output_file.write("DATASET POLYDATA\n")
    output_file.write("POINTS " + str(points) + " float\n")
    for i in range(0, n):
        output_file.write(
            str(p1[0] + radius * np.cos(i * alpha)) + " " + str(p1[1] + radius * np.sin(i * alpha)) + " " +
            str(p1[2]) + "\n")
    for i in range(0, n):
        output_file.write(
            str(p1[0] + radius * np.cos(i * alpha)) + " " + str(p1[1] + radius * np.sin(i * alpha)) + " " +
            str(p1[2] + length * axisZ) + "\n")
    output_file.write(
        "POLYGONS " + str(n + 2) + " " + str(5 * n + (n + 1) * 2) + "\n")

    # writing bottom "circle"
    output_file.write(str(n) + " ")
    for i in range(0, n - 1):
        output_file.write(str(i) + " ")
    output_file.write(str(n - 1) + "\n")

    # writing top "circle"
    output_file.write(str(n) + " ")
    for i in range(0, n - 1):
        output_file.write(str(i + n) + " ")
    output_file.write(str(2 * n - 1) + "\n")

    # writing sides - rectangles
    for i in range(0, n - 1):
        output_file.write("4 " + str(i) + " " + str(
            i + 1) + " " + str(i + n + 1) + " " + str(i + n) + "\n")
    output_file.write("4 " + str(n - 1) + " " + str(
        0) + " " + str(n) + " " + str(2 * n - 1) + "\n")

    output_file.close()
    return 0


def output_vtk_lines(lines, out_file):
    """
    Outputs the VTK files for visualisation of lines in e.g. Paraview.

    Parameters
    ----------
    lines : array_like :obj:`float`
          lines is a list of pairs of points p1, p2
          each pair represents a line segment to output to vtk
          each line in lines contains 6 floats: p1x, p1y, p1z, p2x, p2y, p2z
    out_file : :obj:`str`
          filename for the output

    """

    n_lines = len(lines)

    output_file = open(out_file, "w")
    output_file.write("# vtk DataFile Version 3.0\n")
    output_file.write("Data\n")
    output_file.write("ASCII\n")
    output_file.write("DATASET POLYDATA\n")
    output_file.write("POINTS " + str(2 * n_lines) + " float\n")
    for i in range(0, n_lines):
        one_line = lines[i]
        output_file.write(str(one_line[0]) + " " + str(
            one_line[1]) + " " + str(one_line[2]) + "\n")
        output_file.write(str(one_line[3]) + " " + str(
            one_line[4]) + " " + str(one_line[5]) + "\n")
    output_file.write("LINES " + str(n_lines) + " " + str(3 * n_lines) + "\n")
    for i in range(0, n_lines):
        output_file.write(
            str(2) + " " + str(2 * i) + " " + str(2 * i + 1) + "\n")

    output_file.close()
    return 0


def output_vtk_pore(
        axis, length, outer_rad_left, outer_rad_right, pos, rad_left, rad_right,
        smoothing_radius, m, out_file):
    """
    Outputs the VTK files for visualisation of a pore in e.g. Paraview.

    Parameters
    ----------
    axis : (3,) array_like of :obj:`float`
          The axis
    length : :obj:`float`
          length of pore
    outer_rad_left : :obj:`float`
          outer left radius of pore
    outer_rad_right : :obj:`float`
          outer right radius of pore
    rad_left : :obj:`float`
          inner left radius of pore
    rad_right : :obj:`float`
          inner right radius of pore
    smoothing_radius : :obj:`float`
          smoothing radius for surface connecting outer and inner radii of the pore
    pos : (3,) array_like of :obj:`float`
          Position of the center of the pore
    m : :obj:`int`
          number of discretization sections
    out_file : :obj:`str`
          filename for the output

    """
    # length is the length of the pore without the smoothing part

    # for now, only axis=(1,0,0) is supported
    # should implement rotation
    # m is sufficient to be 10

    if ".vtk" not in out_file:
        print(
            "output_vtk_pore warning: A file with vtk format will be written without .vtk extension.")

    # n must be even therefore:
    n = 2 * m

    # setting points on perimeter
    alpha = 2 * np.pi / n
    beta = 2 * np.pi / n
    number_of_points = 2 * n * (n / 2 + 1)

    output_file = open(out_file, "w")
    output_file.write("# vtk DataFile Version 3.0\n")
    output_file.write("Data\n")
    output_file.write("ASCII\n")
    output_file.write("DATASET POLYDATA\n")
    output_file.write("POINTS " + str(number_of_points) + " float\n")

    # shift center to the left half torus
    p1 = pos - length / 2 * np.array(axis)

    # points on the left half torus
    for j in range(0, n / 2 + 1):
        for i in range(0, n):
            output_file.write(str(p1[0] - np.sin(j * beta)) + " " +
                              str(p1[1] + (rad_left + smoothing_radius - np.cos(j * beta)) * np.cos(i * alpha)) + " " +
                              str(p1[2] + (rad_left + smoothing_radius - np.cos(j * beta)) * np.sin(i * alpha)) + "\n")
    n_points_left = n * (n / 2 + 1)

    # shift center to the right half torus
    p1 = pos + length / 2 * np.array(axis)

    # points on the right half torus
    for j in range(0, n / 2 + 1):
        for i in range(0, n):
            output_file.write(str(p1[0] + np.sin(j * beta)) + " " +
                              str(p1[1] + (rad_right + smoothing_radius - np.cos(j * beta)) * np.cos(i * alpha)) + " " +
                              str(p1[2] + (rad_right + smoothing_radius - np.cos(j * beta)) * np.sin(i * alpha)) + "\n")

    number_of_rectangles = n * n + 2 * n
    output_file.write("POLYGONS " + str(number_of_rectangles)
                      + " " + str(5 * number_of_rectangles) + "\n")

    # writing inner side rectangles
    for i in range(0, n - 1):
        output_file.write("4 " + str(i) + " " +
                          str(i + 1) + " " +
                          str(i + n_points_left + 1) + " " +
                          str(i + n_points_left) + "\n")
    output_file.write("4 " + str(n - 1) + " " +
                      str(0) + " " +
                      str(n_points_left) + " " +
                      str(n_points_left + n - 1) + "\n")

    # writing outer side rectangles
    for i in range(0, n - 1):
        output_file.write("4 " + str(n_points_left - n + i) + " " +
                          str(n_points_left - n + i + 1) + " " +
                          str(n_points_left - n + i + n_points_left + 1) + " " +
                          str(n_points_left - n + i + n_points_left) + "\n")
    output_file.write("4 " + str(n_points_left - n + n - 1) + " " +
                      str(n_points_left - n) + " " +
                      str(n_points_left - n + n_points_left) + " " +
                      str(n_points_left - n + n_points_left + n - 1) + "\n")

    # writing rectangles on the left half of the torus
    for j in range(0, n / 2):
        for i in range(0, n - 1):
            output_file.write("4 " + str(n * j + i) + " " +
                              str(n * j + i + 1) + " " +
                              str(n * j + i + n + 1) + " " +
                              str(n * j + i + n) + "\n")
        output_file.write("4 " + str(n * j + n - 1) + " " +
                          str(n * j) + " " +
                          str(n * j + n) + " " +
                          str(n * j + 2 * n - 1) + "\n")

    # writing rectangles on the right half of the torus
    for j in range(0, n / 2):
        for i in range(0, n - 1):
            output_file.write("4 " + str(n_points_left + n * j + i) + " " +
                              str(n_points_left + n * j + i + 1) + " " +
                              str(n_points_left + n * j + i + n + 1) + " " +
                              str(n_points_left + n * j + i + n) + "\n")
        output_file.write("4 " + str(n_points_left + n * j + n - 1) + " " +
                          str(n_points_left + n * j) + " " +
                          str(n_points_left + n * j + n) + " " +
                          str(n_points_left + n * j + 2 * n - 1) + "\n")

    output_file.close()
    return 0
