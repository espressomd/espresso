import numpy as np
import math

small_epsilon = 0.000000001
large_number = 10000000.0
pi = 3.1415926535897931
output_precision = 14


def custom_str(realn):
    return str('{:.{prec}f}'.format(realn, prec = output_precision))


def get_triangle_normal(a, b, c):
    """
    Returns the normal vector of a triangle given by points a,b,c.

    """
    n = [0.0, 0.0, 0.0]
    n[0] = (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1])
    n[1] = (b[2] - a[2]) * (c[0] - a[0]) - (b[0] - a[0]) * (c[2] - a[2])
    n[2] = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    return np.array(n)


def norm(vect):
    """
    Returns the norm of a vector.

    """
    v = np.array(vect)
    return np.sqrt(np.dot(v,v))


def vec_distance(a, b):
    """
    Returns the length of vector between points a and b.

    """
    return norm(np.array(a) - np.array(b))


def area_triangle(a, b, c):
    """
    Returns the area of a triangle given by points a,b,c.

    """
    n = get_triangle_normal(a, b, c)
    area = 0.5 * norm(n)
    return area


def angle_btw_triangles(P1, P2, P3, P4):
    """
    Returns the size of an angle between triangles given by points P2, P1, P3 and P2, P3, P4.

    """
    n1 = get_triangle_normal(P2, P1, P3)
    n2 = get_triangle_normal(P2, P3, P4)
    tmp11 = np.dot(n1, n2)
    tmp11 = tmp11 * abs(tmp11)

    tmp22 = np.dot(n1, n1)
    tmp33 = np.dot(n2, n2)
    tmp11 /= (tmp22 * tmp33)

    if tmp11 > 0:
        tmp11 = np.sqrt(tmp11)
    else:
        tmp11 = - np.sqrt(- tmp11)

    if tmp11 >= 1.0:
        tmp11 = 0.0
    elif tmp11 <= -1.:
        tmp11 = pi

    phi = pi - math.acos(tmp11)

    if (np.dot(n1, np.array(P4)) - np.dot(n1, np.array(P1))) < 0:
        phi = 2.0 * pi - phi
    return phi


def discard_epsilon(x):
    """
    Returns zero if the argument is too small.

    """
    if (x > -small_epsilon and x < small_epsilon):
        res = 0.0
    else:
        res = x
    return res


def oif_neo_hookean_nonlin(lambd):
    """
    Defines NeoHookean nonlinearity.

    """
    # Defined by (19) from Dupin2007
    res = (pow(lambd, 0.5) + pow(lambd, -2.5)) / (lambd + pow(lambd, -3.))
    return res


def oif_calc_stretching_force(ks, pA, pB, dist0, dist):
    """
    Calculates nonlinear stretching forces between two points on an edge.

    """
    # this has to correspond to the calculation in oif_local_forces.hpp: calc_oif_local
    # as of now, corresponds to git commit f156f9b44dcfd3cef9dd5537a1adfc903ac4772a
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

    """
    # this has to correspond to the calculation in oif_local_forces.hpp: calc_oif_local
    # as of now, corresponds to git commit f156f9b44dcfd3cef9dd5537a1adfc903ac4772a
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

    """
    # this has to correspond to the calculation in oif_local_forces.hpp: calc_oif_local
    # except for division by 3 - each triangle enters this calculation once, while each triangle enters the
    # calc_oif_local three times
    # as of now, corresponds to git commit f156f9b44dcfd3cef9dd5537a1adfc903ac4772a

    centroid = np.array((pA + pB + pC) / 3.0)
    delta_area = A - A0
    ta = centroid - pA
    ta_norm = norm(ta)
    tb = centroid - pB
    tb_norm = norm(tb)
    tc = centroid - pC
    tc_norm = norm(tc)
    common_factor = kal * delta_area / (ta_norm * ta_norm + tb_norm * tb_norm + tc_norm * tc_norm)

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

    """
    # this has to correspond to the calculation in oif_global_forces.hpp: add_oif_global_forces
    # as of now, corresponds to git commit f156f9b44dcfd3cef9dd5537a1adfc903ac4772a
    centroid = np.array((pA + pB + pC) / 3.0)
    delta = Ag - Ag0
    ta = centroid - pA
    ta_norm = norm(ta)
    tb = centroid - pB
    tb_norm = norm(tb)
    tc = centroid - pC
    tc_norm = norm(tc)
    A = area_triangle(pA, pB, pC)
    common_factor = kag * A * delta / (ta_norm * ta_norm + tb_norm * tb_norm + tc_norm * tc_norm)

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

    """
    # this has to correspond to the calculation in oif_global_forces.hpp: add_oif_global_forces
    # as of now, corresponds to git commit f156f9b44dcfd3cef9dd5537a1adfc903ac4772a
    n = get_triangle_normal(pA, pB, pC)
    dn = norm(n)
    vv = (V - V0) / V0
    A = area_triangle(pA, pB, pC)
    f = kv * vv * A * np.array(n) / (dn * 3.0)
    return f


def output_vtk_rhomboid(corner, a, b, c, out_file):
    """
    Outputs the VTK files for visualisation of a rhomboid in e.g. Paraview.

    """
    if ".vtk" not in out_file:
        print("output_vtk_rhomboid warning: A file with vtk format will be written without .vtk extension.")
    output_file = open(out_file, "w")
    output_file.write("# vtk DataFile Version 3.0\n")
    output_file.write("Data\n")
    output_file.write("ASCII\n")
    output_file.write("DATASET POLYDATA\n")
    output_file.write("POINTS 8 float\n")

    output_file.write(str(corner[0]) + " " + str(corner[1]) + " " + str(corner[2]) + "\n")
    output_file.write(str(corner[0] + a[0]) + " " + str(corner[1] + a[1]) + " " + str(corner[2] + a[2]) + "\n")
    output_file.write(str(corner[0] + a[0] + b[0]) + " " + str(corner[1] + a[1] + b[1]) + " " +
                      str(corner[2] + a[2] + b[2]) + "\n")
    output_file.write(str(corner[0] + b[0]) + " " + str(corner[1] + b[1]) + " " + str(corner[2] + b[2]) + "\n")

    output_file.write(str(corner[0] + c[0]) + " " + str(corner[1] + c[1]) + " " + str(corner[2] + c[2]) + "\n")
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


def output_vtk_cylinder(center, axis, length, radius, n, out_file):
    """
    Outputs the VTK files for visualisation of a cylinder in e.g. Paraview.

    """
    # L is the half height of the cylinder
    # only vertical cylinders are supported for now, i.e. with normal (0.0, 0.0, 1.0)

    if ".vtk" not in out_file:
        print("output_vtk_cylinder warning: A file with vtk format will be written without .vtk extension.")
    check_axis = True
    if axis[0]!=0.0:
        check_axis = False
    if axis[1]!=0.0:
        check_axis = False
    if axis[2]==0.0:
        check_axis = False
    if check_axis is False:
        print("output_vtk_cylinder: Output for this type of cylinder is not supported yet.")
        return
    axisZ = 1.0

    # setting points on perimeter
    alpha = 2 * pi / n
    points = 2 * n

    # shift center to the bottom circle
    p1 = center - length * np.array(axis)

    output_file = open(out_file, "w")
    output_file.write("# vtk DataFile Version 3.0\n")
    output_file.write("Data\n")
    output_file.write("ASCII\n")
    output_file.write("DATASET POLYDATA\n")
    output_file.write("POINTS " + str(points) + " float\n")
    for i in range(0,n):
        output_file.write(str(p1[0] + radius*np.cos(i*alpha)) + " " + str(p1[1] + radius*np.sin(i*alpha)) + " " +
                          str(p1[2]) + "\n")
    for i in range(0,n):
        output_file.write(str(p1[0] + radius*np.cos(i*alpha)) + " " + str(p1[1] + radius*np.sin(i*alpha)) + " " +
                          str(p1[2]+2*length*axisZ) + "\n")
    output_file.write("POLYGONS " + str(n+2) + " " + str(5*n+(n+1)*2) + "\n")

    # writing bottom "circle"
    output_file.write(str(n) + " ")
    for i in range(0, n - 1):
        output_file.write(str(i) + " ")
    output_file.write(str(n-1) + "\n")

    # writing top "circle"
    output_file.write(str(n) + " ")
    for i in range(0, n - 1):
        output_file.write(str(i+n) + " ")
    output_file.write(str(2*n - 1) + "\n")

    # writing sides - rectangles
    for i in range(0, n - 1):
        output_file.write("4 " + str(i) + " " + str(i+1) + " " + str(i+n+1) + " " + str(i+n) + "\n")
    output_file.write("4 " + str(n-1) + " " + str(0) + " " + str(n) + " " + str(2*n-1) + "\n")

    output_file.close()
    return 0


def output_vtk_lines(lines, out_file):
    """
    Outputs the VTK files for visualisation of lines in e.g. Paraview.

    """
    # lines is a list of pairs of points p1, p2
    # each pair represents a line segment to output to vtk
    # each line in lines contains 6 floats: p1x, p1y, p1z, p2x, p2y, p2z

    if ".vtk" not in out_file:
        print("output_vtk_lines warning: A file with vtk format will be written without .vtk extension.")

    n_lines = len(lines)

    output_file = open(out_file, "w")
    output_file.write("# vtk DataFile Version 3.0\n")
    output_file.write("Data\n")
    output_file.write("ASCII\n")
    output_file.write("DATASET POLYDATA\n")
    output_file.write("POINTS " + str(2*n_lines) + " float\n")
    for i in range(0,n_lines):
        one_line = lines[i]
        output_file.write(str(one_line[0]) + " " + str(one_line[1]) + " " + str(one_line[2]) + "\n")
        output_file.write(str(one_line[3]) + " " + str(one_line[4]) + " " + str(one_line[5]) + "\n")
    output_file.write("LINES " + str(n_lines) + " " + str(3 *n_lines) + "\n")
    for i in range(0,n_lines):
        output_file.write(str(2) + " " + str(2*i) + " " + str(2*i+1) + "\n")

    output_file.close()
    return 0


def get_lb_interpolated_velocity(position, lbf, system, fluid_agrid):
    # position is a vector [x,y,z], in which we want the fluid velocity

    # since the fluid node [0,0,0] is at the position [0.5, 0.5, 0.5], we need to shift, fold and rescale
    shifted_position = position - np.array(0.5*fluid_agrid, 0.5*fluid_agrid, 0.5*fluid_agrid)
    shifted_position = np.mod(shifted_position, system.box_l)/fluid_agrid

    # these six variables are not positions, but lattice indices
    x_lower = math.floor(shifted_position[0])
    x_upper = math.ceil(shifted_position[0])
    y_lower = math.floor(shifted_position[1])
    y_upper = math.ceil(shifted_position[1])
    z_lower = math.floor(shifted_position[2])
    z_upper = math.ceil(shifted_position[2])

    vel000 = np.array(lbf[x_lower, y_lower, z_lower].velocity)
    vel100 = np.array(lbf[x_upper, y_lower, z_lower].velocity)
    vel010 = np.array(lbf[x_lower, y_upper, z_lower].velocity)
    vel001 = np.array(lbf[x_lower, y_lower, z_upper].velocity)
    vel110 = np.array(lbf[x_upper, y_upper, z_lower].velocity)
    vel101 = np.array(lbf[x_upper, y_lower, z_upper].velocity)
    vel011 = np.array(lbf[x_lower, y_upper, z_upper].velocity)
    vel111 = np.array(lbf[x_upper, y_upper, z_upper].velocity)

    x = shifted_position[0] - x_lower
    y = shifted_position[1] - y_lower
    z = shifted_position[2] - z_lower

    interpolated_velocity = vel000 * (1 - x) * (1 - y) * (1 - z) \
        + vel001 * (1 - x) * (1 - y) * z \
        + vel010 * (1 - x) * y * (1 - z) \
        + vel100 * x * (1 - y) * (1 - z) \
        + vel110 * x * y * (1 - z) \
        + vel101 * x * (1 - y) * z \
        + vel011 * (1 - x) * y * z \
        + vel111 * x * y * z

    return interpolated_velocity
