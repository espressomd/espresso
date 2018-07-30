import espressomd
import numpy as np
from .oif_utils import *
from espressomd.interactions import OifLocalForces
from espressomd.interactions import OifGlobalForces
from espressomd.interactions import OifOutDirection


class FixedPoint(object):
    """
    Represents mesh points, not connected to any ESPResSo particle.

    """
    def __init__(self, pos, id):
        if not isinstance(id, int):
            raise TypeError("Id must be integer.")
        if not ((len(pos) == 3) and isinstance(pos[0],float) and isinstance(pos[1],float) and isinstance(pos[2],float)):
            raise TypeError("Pos must be a list of three floats.") 

        self.x = pos[0]
        self.y = pos[1]
        self.z = pos[2]
        self.id = id


    def get_pos(self):
        return [self.x, self.y, self.z]

    def get_id(self):
        return self.id


class PartPoint(object):
    """
    Represents mesh points, connected to ESPResSo particle.

    """
    def __init__(self, part, id, part_id):  # part is physical ESPResSo particle corresponding to that particular point
        if not (isinstance(part, espressomd.particle_data.ParticleHandle) and isinstance(id,int) and isinstance(part_id,int)):
            raise TypeError("Arguments to PartPoint are incorrect.")
        self.part = part
        self.part_id = part_id  # because in adding bonds to the particles in OifCell
        # one needs to know the global id of the particle.
        self.id = id

    def get_pos(self):
        return self.part.pos

    def get_vel(self):
        return self.part.v

    def get_mass(self):
        return self.part.mass

    def get_type(self):
        return self.part.type

    def get_force(self):
        return self.part.f

    def set_pos(self,pos):
        self.part.pos = pos

    def set_vel(self, vel):
        self.part.v = vel

    def set_force(self, force):
        self.part.ext_force = force

    def kill_motion(self):
        self.part.fix = [1, 1, 1]
        
    def unkill_motion(self):
        self.part.unfix()


class Edge(object):
    """
    Represents edges in a mesh.

    """
    def __init__(self, A, B):
        if not (isinstance(A,PartPoint) or (isinstance(A,FixedPoint))) and (isinstance(B,PartPoint) or (isinstance(B,FixedPoint))):
            TypeError("Arguments to Edge must be FixedPoint or PartPoint.")
        self.A = A
        self.B = B

    def length(self):
        return vec_distance(self.A.get_pos(), self.B.get_pos())


class Triangle(object):
    """
    Represents triangles in a mesh.

    """
    def __init__(self, A, B, C):
        if not (isinstance(A,PartPoint) or (isinstance(A,FixedPoint))) and (isinstance(B,PartPoint) or (isinstance(B,FixedPoint))) and (isinstance(C,PartPoint) or (isinstance(C,FixedPoint))):
            TypeError("Arguments to Triangle must be FixedPoint or PartPoint.")
        self.A = A
        self.B = B
        self.C = C

    def area(self):
        area = area_triangle(self.A.get_pos(), self.B.get_pos(), self.C.get_pos())
        return area


class Angle(object):
    """
    Represents angles in a mesh.

    """
    def __init__(self, A, B, C, D):
        if not (isinstance(A, PartPoint) or (isinstance(A, FixedPoint))) \
                and (isinstance(B, PartPoint) or (isinstance(B, FixedPoint))) \
                and (isinstance(C, PartPoint) or (isinstance(C, FixedPoint))) \
                and (isinstance(D, PartPoint) or (isinstance(D, FixedPoint))):
            TypeError("Arguments to Angle must be FixedPoint or PartPoint.")
        self.A = A
        self.B = B
        self.C = C
        self.D = D

    def size(self):
        angle_size = angle_btw_triangles(self.A.get_pos(), self.B.get_pos(), self.C.get_pos(), self.D.get_pos())
        return angle_size


class ThreeNeighbors(object):
    """
    Represents three best spatially distributed neighbors of a point in a mesh.

    """
    def __init__(self, A, B, C):
        if not (isinstance(A, PartPoint) or (isinstance(A, FixedPoint))) \
                and (isinstance(B, PartPoint) or (isinstance(B, FixedPoint))) \
                and (isinstance(C, PartPoint) or (isinstance(C, FixedPoint))):
            TypeError("Arguments to ThreeNeighbors must be FixedPoint or PartPoint.")
        self.A = A
        self.B = B
        self.C = C
   

    def outer_normal(self):
        outer_normal = get_triangle_normal(self.A.get_pos(), self.B.get_pos(), self.C.get_pos())
        return outer_normal


class Mesh(object):
    """
    Represents a triangular mesh.

    """
    def __init__(self, nodes_file=None, triangles_file=None, system=None, resize=(1.0, 1.0, 1.0),
                 particle_type=-1, particle_mass=1.0, normal=False, check_orientation=True):
        if (system is None) or (not isinstance(system,espressomd.System)):
            raise Exception("Mesh: No system provided or wrong type given. Quitting.")

        self.system = system
        self.normal = normal
        self.nodes_file = nodes_file
        self.triangles_file = triangles_file

        self.points = []
        self.edges = []
        self.triangles = []
        self.angles = []
        self.neighbors = []
        self.ids_extremal_points = [0, 0, 0, 0, 0, 0, 0]

        if not ((nodes_file is None) or (triangles_file is None)):
            if not (isinstance(nodes_file,str) and isinstance(triangles_file,str)):
                raise TypeError("Mesh: Filenames must be strings.")
            if not ((len(resize) == 3) and isinstance(resize[0],float) and isinstance(resize[1],float) and isinstance(resize[2],float)):
                raise TypeError("Mesh: Pos must be a list of three floats.") 
            if not isinstance(particle_type,int):
                raise TypeError("Mesh: particle_type must be integer.")
            if not isinstance(particle_mass,float):
                raise TypeError("Mesh: particle_mass must be float.")
            if not isinstance(normal,bool):
                raise TypeError("Mesh: normal must be bool.") 
            if not isinstance(check_orientation,bool):
                raise TypeError("Mesh: check_orientation must be bool.") 
            # reading the mesh point positions from file
            in_file = open(nodes_file, "r")
            nodes_coord = in_file.read().split("\n")
            in_file.close()
            # removes a blank line at the end of the file if there is any:
            nodes_coord = filter(None, nodes_coord) # here we have list of lines with triplets of strings
            for line in nodes_coord:  # extracts coordinates from the string line
                line = np.array([float(x) for x in line.split()])
                coords = np.array(resize) * line
                tmp_fixed_point = FixedPoint(coords, len(self.points))
                self.points.append(tmp_fixed_point)
    
            # searching for extremal points IDs
            x_min = large_number
            x_max = -large_number
            y_min = large_number
            y_max = -large_number
            z_min = large_number
            z_max = -large_number
            for tmp_fixed_point in self.points:
                coords = tmp_fixed_point.get_pos()
                if coords[0] < x_min:
                    x_min = coords[0]
                    self.ids_extremal_points[0] = tmp_fixed_point.get_id()
                if coords[0] > x_max:
                    x_max = coords[0]
                    self.ids_extremal_points[1] = tmp_fixed_point.get_id()
                if coords[1] < y_min:
                    y_min = coords[1]
                    self.ids_extremal_points[2] = tmp_fixed_point.get_id()
                if coords[1] > y_max:
                    y_max = coords[1]
                    self.ids_extremal_points[3] = tmp_fixed_point.get_id()
                if coords[2] < z_min:
                    z_min = coords[2]
                    self.ids_extremal_points[4] = tmp_fixed_point.get_id()
                if coords[2] > z_max:
                    z_max = coords[2]
                    self.ids_extremal_points[5] = tmp_fixed_point.get_id()
    
            # reading the triangle incidences from file
            in_file = open(triangles_file, "r")
            triangles_incid = in_file.read().split("\n")
            in_file.close()
            # removes a blank line at the end of the file if there is any:
            triangles_incid = filter(None, triangles_incid)
            for line in triangles_incid:  # extracts incidences from the string line
                incid = np.array([int(x) for x in line.split()])
                tmp_triangle = Triangle(self.points[incid[0]], self.points[incid[1]], self.points[incid[2]])
                self.triangles.append(tmp_triangle)
    
            if check_orientation is True:
                # check whether all triangles in file had the same orientation; if not, correct the orientation
                self.check_orientation()
    
            # creating list of edge incidences from triangle incidences
            # using temporary list of edge incidences
            tmp_edge_incidences = []
            for triangle in self.triangles:
                pa = triangle.A.id
                pb = triangle.B.id
                pc = triangle.C.id
                if ([pa, pb] not in tmp_edge_incidences) and ([pb, pa] not in tmp_edge_incidences):
                    tmp_edge_incidences.append([pa, pb])
                if ([pb, pc] not in tmp_edge_incidences) and ([pc, pb] not in tmp_edge_incidences):
                    tmp_edge_incidences.append([pb, pc])
                if ([pa, pc] not in tmp_edge_incidences) and ([pc, pa] not in tmp_edge_incidences):
                    tmp_edge_incidences.append([pa, pc])
            for tmp_incid in tmp_edge_incidences:
                tmp_edge = Edge(self.points[tmp_incid[0]], self.points[tmp_incid[1]])
                self.edges.append(tmp_edge)
    
            # creating list angles (former bending incidences) from triangle incidences
            for edge in self.edges:
                pa = edge.A.id
                pb = edge.B.id
                pc = -1
                pd = -1
                detected = 0
                # detected = number of detected triangles with current edge common
                # Algorithm is as follows: we run over all triangles and check
                # whether two vertices are those from current edge. If we find such triangle,
                # we put the ID of the third vertex to pc and we check if the orientation pa, pb, pc is the same as
                # was in the triangle list (meaning, that we found one of the following three triples
                # in the triangle list: pa, pb, pc or pb, pc, pa or pc, pa, pb).
                # If we have the same orientation, we set orient = 1, otherwise orient = -1.
                # Then we go further looking for the second triangle.
                # The second triangle should have the opposite orientation.
                # The normal of the first triangle will be P1P2 x P1P3, of the second triangle will be P2P4 x P2P3
                orient = 0
                for triangle in self.triangles:
                    # Run over all triangles and determine the two triangles with the common current edge
                    if (pa == triangle.A.id) and (pb == triangle.B.id):
                        if detected == 0:
                            # if no triangle with such edge was detected before
                            pc = triangle.C.id
                            detected = 1
                            orient = 1
                        else:
                            # if this is the second triangle with this edge, then also quit the for-loop over triangles 
                            pd = triangle.C.id
                            break
                    if (pa == triangle.B.id) and (pb == triangle.C.id):
                        if detected == 0:
                            pc = triangle.A.id
                            detected = 1
                            orient = 1
                        else:
                            pd = triangle.A.id
                            break
                    if (pa == triangle.C.id) and (pb == triangle.A.id):
                        if detected == 0:
                            pc = triangle.B.id
                            detected = 1
                            orient = 1
                        else:
                            pd = triangle.B.id
                            break
                    if (pa == triangle.B.id) and (pb == triangle.A.id):
                        if detected == 0:
                            pc = triangle.C.id
                            detected = 1
                            orient = -1
                        else:
                            pd = triangle.C.id
                            break
                    if (pa == triangle.C.id) and (pb == triangle.B.id):
                        if detected == 0:
                            pc = triangle.A.id
                            detected = 1
                            orient = -1
                        else:
                            pd = triangle.A.id
                            break
                    if (pa == triangle.A.id) and (pb == triangle.C.id):
                        if detected == 0:
                            pc = triangle.B.id
                            detected = 1
                            orient = -1
                        else:
                            pd = triangle.B.id
                            break
                if orient == 1:
                    tmp = pd
                    pd = pc
                    pc = tmp
                tmp_angle = Angle(self.points[pc], self.points[pa], self.points[pb], self.points[pd])
                self.angles.append(tmp_angle)
    
            # creating list of three neighbors for membrane collision
            if normal is True:
                for point in self.points:
                    tmp_neighbors = []
                    # cycle through edges and select those that contain point
                    for edge in self.edges:
                        # take an edge and copy the nodes of the edge to pa, pb
                        if edge.A.id == point.id:
                            tmp_neighbors.append(edge.B)
                        if edge.B.id == point.id:
                            tmp_neighbors.append(edge.A)
                    # create vectors to all neighbors and normalize them
                    tmp_vectors_to_neighbors = []
                    p_coords = np.array(point.get_pos())
                    for neighbor in tmp_neighbors:
                        tmp_vector = neighbor.get_pos() - p_coords
                        tmp_length = norm(tmp_vector)
                        if tmp_length < small_epsilon:
                            raise Exception("Mesh: Degenerate edge. Quitting.")
                        tmp_vector /= tmp_length
                        tmp_vectors_to_neighbors.append(tmp_vector)
                    # check all triplets of neighbors and select the one that is best spatially distributed
                    # by adding the corresponding three normalized vectors
                    # and selecting the one with smallest resultant vector
                    n_neighbors = len(tmp_neighbors)
                    min_length = large_number
                    best_neighbors = [tmp_neighbors[0], tmp_neighbors[1], tmp_neighbors[2]]
                    for i in range(0,n_neighbors):
                        for j in range(i+1,n_neighbors):
                            for k in range(j+1,n_neighbors):
                                tmp_result_vector = tmp_vectors_to_neighbors[i] + tmp_vectors_to_neighbors[j] + \
                                                    tmp_vectors_to_neighbors[k]
                                tmp_result_vector_length = norm(tmp_result_vector)
                                if tmp_result_vector_length < min_length:
                                    min_length = tmp_result_vector_length
                                    best_neighbors = [tmp_neighbors[i], tmp_neighbors[j], tmp_neighbors[k]]
                    # find one triangle that contains this point and compute its normal vector
                    for triangle in self.triangles:
                        if triangle.A.id == point.id or triangle.B.id == point.id or triangle.C.id == point.id:
                            tmp_normal_triangle = get_triangle_normal(triangle.A.get_pos(), triangle.B.get_pos(),
                                                                      triangle.C.get_pos())
                            break
                    # properly orient selected neighbors and save them to the list of neighbors
                    tmp_normal_neighbors = get_triangle_normal(best_neighbors[0].get_pos(), best_neighbors[1].get_pos(),
                                                               best_neighbors[2].get_pos())
                    tmp_length_normal_triangle = norm(tmp_normal_triangle)
                    tmp_length_normal_neighbors = norm(tmp_normal_neighbors)
                    tmp_product = np.dot(tmp_normal_triangle, tmp_normal_neighbors) / \
                                  (tmp_length_normal_triangle * tmp_length_normal_neighbors)
                    tmp_angle = np.arccos(tmp_product)
                    if tmp_angle > np.pi/2.0:
                        selected_neighbors = ThreeNeighbors(best_neighbors[0], best_neighbors[1], best_neighbors[2])
                    else:
                        selected_neighbors = ThreeNeighbors(best_neighbors[0], best_neighbors[2], best_neighbors[1])
                    self.neighbors.append(selected_neighbors)
            else:
                for point in self.points:
                    selected_neighbors = ThreeNeighbors(point, point, point)
                    self.neighbors.append(selected_neighbors)

    def copy(self, origin=None, particle_type=-1, particle_mass=1.0, rotate=None):
        mesh = Mesh(system=self.system)
        mesh.ids_extremal_points = self.ids_extremal_points
        rotation = np.array([[1.0, 0.0, 0.0],
                            [0.0, 1.0, 0.0],
                            [0.0, 0.0, 1.0]])

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
            tmp_pos = point.get_pos()
            tmp_rotate_pos = np.array(point.get_pos())
            # rotation of nodes
            if rotate is not None:
                tmp_pos = rotation.dot(tmp_rotate_pos)
                tmp_pos = [discard_epsilon(tmp_pos[0]), discard_epsilon(tmp_pos[1]), discard_epsilon(tmp_pos[2])]
            if origin is not None:
                tmp_pos += np.array(origin)
            new_part_id = len(self.system.part)  # to remember the global id of the ESPResSo particle
            self.system.part.add(pos=tmp_pos, type=particle_type, mass=particle_mass, mol_id=particle_type)
            new_part = self.system.part[new_part_id]
            new_part_point = PartPoint(new_part, len(mesh.points), new_part_id)
            mesh.points.append(new_part_point)
        for edge in self.edges:
            new_edge = Edge(mesh.points[edge.A.id], mesh.points[edge.B.id])
            mesh.edges.append(new_edge)
        for triangle in self.triangles:
            new_triangle = Triangle(mesh.points[triangle.A.id], mesh.points[triangle.B.id], mesh.points[triangle.C.id])
            mesh.triangles.append(new_triangle)
        for angle in self.angles:
            new_angle = Angle(mesh.points[angle.A.id], mesh.points[angle.B.id], mesh.points[angle.C.id],
                              mesh.points[angle.D.id])
            mesh.angles.append(new_angle)
        for neighbors in self.neighbors:
            new_neighbors = ThreeNeighbors(mesh.points[neighbors.A.id], mesh.points[neighbors.B.id],
                                           mesh.points[neighbors.C.id])
            mesh.neighbors.append(new_neighbors)
        return mesh

    def check_orientation(self):
        tmp_triangle_list = []
        tmp_triangle_list_ok = []
        t_ok = None
        corrected_triangle = None
        for triangle in self.triangles:
            tmp_triangle_list.append(triangle)

        # move the first triangle to the checked and corrected list
        tmp_triangle_list_ok.append(tmp_triangle_list[0])
        tmp_triangle_list.pop(0)

        while len(tmp_triangle_list) != 0:
            i = 0
            while i < len(tmp_triangle_list):
                tmp_triangle = tmp_triangle_list[i]
                for correct_triangle in tmp_triangle_list_ok:
                    # check if triangles have a common edge, if so, check orientation
                    are_neighbors = True
                    if tmp_triangle.A.id == correct_triangle.A.id:
                        if tmp_triangle.B.id == correct_triangle.B.id:
                            t_ok = False  # this is situation 123 and 124
                            corrected_triangle = Triangle(tmp_triangle.A, tmp_triangle.C, tmp_triangle.B)
                        else:
                            if tmp_triangle.B.id == correct_triangle.C.id:
                                t_ok = True  # this is situation 123 and 142
                            else:
                                if tmp_triangle.C.id == correct_triangle.B.id:
                                    t_ok = True  # this is situation 123 and 134
                                else:
                                    if tmp_triangle.C.id == correct_triangle.C.id:
                                        t_ok = False  # this is situation 123 and 143
                                        corrected_triangle = Triangle(tmp_triangle.A, tmp_triangle.C, tmp_triangle.B)
                                    else:
                                        are_neighbors = False
                    else:
                        if tmp_triangle.A.id == correct_triangle.B.id:
                            if tmp_triangle.B.id == correct_triangle.C.id:
                                t_ok = False  # this is situation 123 and 412
                                corrected_triangle = Triangle(tmp_triangle.A, tmp_triangle.C, tmp_triangle.B)
                            else:
                                if tmp_triangle.B.id == correct_triangle.A.id:
                                    t_ok = True  # this is situation 123 and 214
                                else:
                                    if tmp_triangle.C.id == correct_triangle.C.id:
                                        t_ok = True  # this is situation 123 and 413
                                    else:
                                        if tmp_triangle.C.id == correct_triangle.A.id:
                                            t_ok = False  # this is situation 123 and 314
                                            corrected_triangle = Triangle(tmp_triangle.A, tmp_triangle.C,
                                                                          tmp_triangle.B)
                                        else:
                                            are_neighbors = False
                        else:
                            if tmp_triangle.A.id == correct_triangle.C.id:
                                if tmp_triangle.B.id == correct_triangle.A.id:
                                    t_ok = False  # this is situation 123 and 241
                                    corrected_triangle = Triangle(tmp_triangle.A, tmp_triangle.C, tmp_triangle.B)
                                else:
                                    if tmp_triangle.B.id == correct_triangle.B.id:
                                        t_ok = True  # this is situation 123 and 421
                                    else:
                                        if tmp_triangle.C.id == correct_triangle.A.id:
                                            t_ok = True  # this is situation 123 and 341
                                        else:
                                            if tmp_triangle.C.id == correct_triangle.B.id:
                                                t_ok = False  # this is situation 123 and 431
                                                corrected_triangle = Triangle(tmp_triangle.A, tmp_triangle.C,
                                                                              tmp_triangle.B)
                                            else:
                                                are_neighbors = False
                            else:
                                if tmp_triangle.B.id == correct_triangle.A.id:
                                    if tmp_triangle.C.id == correct_triangle.B.id:
                                        t_ok = False  # this is situation 123 and 234
                                        corrected_triangle = Triangle(tmp_triangle.A, tmp_triangle.C, tmp_triangle.B)
                                    else:
                                        if tmp_triangle.C.id == correct_triangle.C.id:
                                            t_ok = True  # this is situation 123 and 243
                                        else:
                                            are_neighbors = False
                                else:
                                    if tmp_triangle.B.id == correct_triangle.B.id:
                                        if tmp_triangle.C.id == correct_triangle.C.id:
                                            t_ok = False  # this is situation 123 and 423
                                            corrected_triangle = Triangle(tmp_triangle.A, tmp_triangle.C,
                                                                          tmp_triangle.B)
                                        else:
                                            if tmp_triangle.C.id == correct_triangle.A.id:
                                                t_ok = True  # this is situation 123 and 324
                                            else:
                                                are_neighbors = False
                                    else:
                                        if tmp_triangle.B.id == correct_triangle.C.id:
                                            if tmp_triangle.C.id == correct_triangle.A.id:
                                                t_ok = False  # this is situation 123 and 342
                                                corrected_triangle = Triangle(tmp_triangle.A, tmp_triangle.C,
                                                                              tmp_triangle.B)
                                            else:
                                                if tmp_triangle.C.id == correct_triangle.B.id:
                                                    t_ok = True  # this is situation 123 and 432
                                                else:
                                                    are_neighbors = False
                                        else:
                                            are_neighbors = False
                    if are_neighbors:
                        # move the tmp_triangle to the checked and corrected list
                        if t_ok:
                            tmp_triangle_list_ok.append(tmp_triangle)
                        else:
                            tmp_triangle_list_ok.append(corrected_triangle)
                        tmp_triangle_list.pop(i)
                        break
                i += 1
        # replace triangles with checked triangles
        i = 0
        for tmp_triangle in tmp_triangle_list_ok:
            self.triangles[i] = Triangle(tmp_triangle.A, tmp_triangle.C, tmp_triangle.B)
            i += 1
        # all triangles now have the same orientation, check if it is correct
        tmp_volume = self.volume()
        if tmp_volume < 0:
            # opposite orientation, flip all triangles
            i = 0
            for tmp_triangle in self.triangles:
                    self.triangles[i] = Triangle(tmp_triangle.A, tmp_triangle.C, tmp_triangle.B)
                    i += 1
        return 0

    def surface(self):
        surface = 0.0
        for triangle in self.triangles:
            surface += triangle.area()
        return surface

    def volume(self):
        volume = 0.0
        for triangle in self.triangles:
            tmp_normal = get_triangle_normal(triangle.A.get_pos(), triangle.B.get_pos(), triangle.C.get_pos())
            tmp_normal_length = norm(tmp_normal)
            tmp_sum_z_coords = 1.0 / 3.0 * (triangle.A.get_pos()[2] + triangle.B.get_pos()[2] + triangle.C.get_pos()[2])
            volume -= triangle.area() * tmp_normal[2] / tmp_normal_length * tmp_sum_z_coords
        return volume

    def get_n_nodes(self):
        return len(self.points)

    def get_n_triangles(self):
        return len(self.triangles)
        
    def get_n_edges(self):
        return len(self.edges)

    def output_mesh_triangles(self, triangles_file=None):
        # this is useful after the mesh correction
        # output of mesh nodes can be done from OifCell (this is because their position may change)
        if triangles_file is None:
            raise Exception("OifMesh: No file_name provided for triangles. Quitting.")
        output_file = open(triangles_file, "w")
        for t in self.triangles:
            output_file.write(str(t.A.id) + " " + str(t.B.id) + " " + str(t.C.id) + "\n")
        output_file.close()
        return 0

    def mirror(self, mirror_x=0, mirror_y=0, mirror_z=0, out_file_name=""):
        if out_file_name == "":
            raise Exception("Cell.Mirror: output meshnodes file for new mesh is missing. Quitting.")
        if (mirror_x!=0 and mirror_x != 1) or (mirror_y!=0 and mirror_y != 1) or (mirror_z!=0 and mirror_z != 1):
            raise Exception("Mesh.Mirror: for mirroring only values 0 or 1 are accepted. 1 indicates that the corresponding coordinate will be flipped.  Exiting.")
        if mirror_x + mirror_y + mirror_z > 1:
            raise Exception("Mesh.Mirror: flipping allowed only for one axis. Exiting.")
        if mirror_x + mirror_y + mirror_z == 1:
            out_file = open(out_file_name, "w")
            for p in self.points:
                coor = p.get_pos()
                if mirror_x == 1:
                    coor[0] *= -1.0
                if mirror_y == 1:
                    coor[1] *= -1.0
                if mirror_z == 1:
                    coor[2] *= -1.0
                out_file.write(custom_str(coor[0]) + " " + custom_str(coor[1]) + " " + custom_str(coor[2]) + "\n")
            out_file.close()
        return 0


class OifCellType(object):  # analogous to oif_template
    """
    Represents a template for creating elastic objects.

    """
    def __init__(self, nodes_file="", triangles_file="", system=None, resize=(1.0, 1.0, 1.0), ks=0.0, kslin=0.0,
                 kb=0.0, kal=0.0, kag=0.0, kv=0.0, kvisc=0.0, normal=False, check_orientation=True):
        if (system is None) or (not isinstance(system,espressomd.System)):
            raise Exception("OifCellType: No system provided or wrong type. Quitting.")
        if (nodes_file == "") or (triangles_file == ""):
            raise Exception("OifCellType: One of nodesfile or trianglesfile is missing. Quitting.")
        if not (isinstance(nodes_file,str) and isinstance(triangles_file,str)):
            raise TypeError("OifCellType: Filenames must be strings.")
        if not ((len(resize) == 3) and isinstance(resize[0],float) and isinstance(resize[1],float) and isinstance(resize[2],float)):
            raise TypeError("OifCellType: Resize must be a list of three floats.") 
        if not (isinstance(ks,float) and isinstance(ks,float) and isinstance(kb,float) and isinstance(kal,float) and isinstance(kag,float) and isinstance(kv,float) and isinstance(kvisc,float)):
            raise TypeError("OifCellType: Elastic parameters must be floats.")
        if not isinstance(normal,bool):
            raise TypeError("OifCellType: normal must be bool.") 
        if not isinstance(check_orientation,bool):
            raise TypeError("OifCellType: check_orientation must be bool.")     
        if (ks != 0.0) and (kslin != 0.0):
            raise Exception("OifCellType: Cannot use linear and nonlinear stretching at the same time. Quitting.")
        self.system = system
        self.mesh = Mesh(nodes_file=nodes_file, triangles_file=triangles_file, system=system, resize=resize,
                         normal=normal, check_orientation=check_orientation)
        self.local_force_interactions = []
        self.resize = resize
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
                r0 = vec_distance(angle.B.get_pos(), angle.C.get_pos())
                phi = angle_btw_triangles(angle.A.get_pos(), angle.B.get_pos(), angle.C.get_pos(), angle.D.get_pos())
                area1 = area_triangle(angle.A.get_pos(), angle.B.get_pos(), angle.C.get_pos())
                area2 = area_triangle(angle.D.get_pos(), angle.B.get_pos(), angle.C.get_pos())
                tmp_local_force_inter = OifLocalForces(r0=r0, ks=ks, kslin=kslin, phi0=phi, kb=kb, A01=area1, A02=area2,
                                                       kal=kal, kvisc=kvisc)
                self.local_force_interactions.append([tmp_local_force_inter, [angle.A, angle.B, angle.C, angle.D]])
                self.system.bonded_inter.add(tmp_local_force_inter)
        if (kag != 0.0) or (kv != 0.0):
            surface = self.mesh.surface()
            volume = self.mesh.volume()
            self.global_force_interaction = OifGlobalForces(A0_g=surface, ka_g=kag, V0=volume, kv=kv)
            self.system.bonded_inter.add(self.global_force_interaction)

    def print_info(self):
        print("\nThe following OifCellType was created: ")
        print("\t nodes_file: " + self.mesh.nodes_file)
        print("\t triangles_file: " + self.mesh.triangles_file)
        print("\t n_nodes: " + str(self.mesh.get_n_nodes()))
        print("\t n_triangles: " + str(self.mesh.get_n_triangles()))
        print("\t n_edges: " + str(self.mesh.get_n_edges()))
        print("\t ks: " + custom_str(self.ks))
        print("\t kslin: " + custom_str(self.kslin))
        print("\t kb: " + custom_str(self.kb))
        print("\t kal: " + custom_str(self.kal))
        print("\t kag: " + custom_str(self.kag))
        print("\t kv: " + custom_str(self.kv))
        print("\t kvisc: " + custom_str(self.kvisc))
        print("\t normal: " + str(self.normal))
        print("\t resize: " + str(self.resize))
        print(" ")


class OifCell(object):
    """
    Represents a concrete elastic object.

    """
    def __init__(self, cell_type=None, origin=None, particle_type=None, particle_mass=1.0, rotate=None):
        if (cell_type is None) or (not isinstance(cell_type,OifCellType)):
            raise Exception("OifCell: No cellType provided or wrong type. Quitting.")
        if (origin is None) or \
        (not ((len(origin) == 3) and isinstance(origin[0],float) and isinstance(origin[1],float) and isinstance(origin[2],float))):
            raise TypeError("Origin must be tuple.")
        if (particle_type is None) or (not isinstance(particle_type,int)):
            raise Exception("OifCell: No particle_type specified or wrong type. Quitting.")
        if not isinstance(particle_mass,float):
            raise Exception("OifCell: particle mass must be float.")
        if (rotate is not None) and not ((len(rotate) == 3) and isinstance(rotate[0],float) and isinstance(rotate[1],float) and isinstance(rotate[2],float)):
            raise TypeError("Rotate must be list of three floats.")

        self.cell_type = cell_type
        self.cell_type.system.max_oif_objects =self.cell_type.system.max_oif_objects+1
        self.mesh = cell_type.mesh.copy(origin=origin, particle_type=particle_type, particle_mass=particle_mass, rotate=rotate)
        self.particle_mass = particle_mass
        self.particle_type = particle_type
        self.origin = origin
        self.rotate = rotate
        for inter in self.cell_type.local_force_interactions:
            esp_inter = inter[0]
            points = inter[1]
            n_points = len(points)
            if n_points == 2:
                p0 = self.mesh.points[points[0].id]  # Getting PartPoints from id's of FixedPoints
                p1 = self.mesh.points[points[1].id]
                p0.part.add_bond((esp_inter, p1.part_id))
            if n_points == 3:
                p0 = self.mesh.points[points[0].id]
                p1 = self.mesh.points[points[1].id]
                p2 = self.mesh.points[points[2].id]
                p0.part.add_bond((esp_inter, p1.part_id, p2.part_id))
            if n_points == 4:
                p0 = self.mesh.points[points[0].id]
                p1 = self.mesh.points[points[1].id]
                p2 = self.mesh.points[points[2].id]
                p3 = self.mesh.points[points[3].id]
                p1.part.add_bond((esp_inter, p0.part_id, p2.part_id, p3.part_id))

        if (self.cell_type.kag!=0.0) or (self.cell_type.kv!=0.0):
            for triangle in self.mesh.triangles:
                triangle.A.part.add_bond((self.cell_type.global_force_interaction, triangle.B.part_id,
                                          triangle.C.part_id))

        # setting the out_direction interaction for membrane collision
        if self.cell_type.mesh.normal is True:
            tmp_out_direction_interaction = OifOutDirection()
            # this interaction could be just one for all objects, but here it is created multiple times
            self.cell_type.system.bonded_inter.add(tmp_out_direction_interaction)
            for p in self.mesh.points:
                p.part.add_bond((tmp_out_direction_interaction, self.mesh.neighbors[p.id].A.part_id,
                                 self.mesh.neighbors[p.id].B.part_id, self.mesh.neighbors[p.id].C.part_id))


    
    def get_origin(self):
        center = np.array([0.0, 0.0, 0.0])
        for p in self.mesh.points:
            center += p.get_pos()
        return center/len(self.mesh.points)

    def set_origin(self, new_origin = (0.0, 0.0, 0.0)):
        old_origin = self.get_origin()
        for p in self.mesh.points:
            new_position = p.get_pos() - old_origin + new_origin
            p.set_pos(new_position)

    def get_approx_origin(self):
        approx_center = np.array([0.0, 0.0, 0.0])
        for id in self.mesh.ids_extremal_points:
            approx_center += self.mesh.points[id].get_pos()
        return approx_center/len(self.mesh.ids_extremal_points)

    def get_origin_folded(self):
        origin = self.get_origin()
        return np.mod(origin, self.cell_type.system.box_l)

    def get_velocity(self):
        velocity = np.array([0.0, 0.0, 0.0])
        for p in self.mesh.points:
            velocity += p.get_vel()
        return velocity/len(self.mesh.points)

    def set_velocity(self, new_velocity = (0.0, 0.0, 0.0)):
        for p in self.mesh.points:
            p.set_vel(new_velocity)

    def pos_bounds(self):
        x_min = large_number
        x_max = -large_number
        y_min = large_number
        y_max = -large_number
        z_min = large_number
        z_max = -large_number
        for p in  self.mesh.points:
            coords = p.get_pos()
            if coords[0] < x_min:
               x_min = coords[0]
            if coords[0] > x_max:
                x_max = coords[0]
            if coords[1] < y_min:
                y_min = coords[1]
            if coords[1] > y_max:
                y_max = coords[1]
            if coords[2] < z_min:
                z_min = coords[2]
            if coords[2] > z_max:
                z_max = coords[2]
        return [x_min, x_max, y_min, y_max, z_min, z_max]

    def surface(self):
        return self.mesh.surface()

    def volume(self):
        return self.mesh.volume()

    def diameter(self):
        max_distance = 0.0
        n_points = len(self.mesh.points)
        for i in range(0, n_points):
            for j in range(i+1, n_points):
                p1 = self.mesh.points[i].get_pos()
                p2 = self.mesh.points[j].get_pos()
                tmp_dist = vec_distance(p1,p2)
                if tmp_dist > max_distance:
                    max_distance = tmp_dist
        return max_distance
        
    def get_n_nodes(self):
        return self.mesh.get_n_nodes()

    def set_force(self, new_force = (0.0, 0.0, 0.0)):
        for p in self.mesh.points:
            p.set_force(new_force)

    # this is not implemented
    # def kill_motion(self):
    #    for p in self.mesh.points:
    #        p.kill_motion()

    # this is not implemented
    # def unkill_motion(self):
    #    for p in self.mesh.points:
    #        p.unkill_motion()

    def output_vtk_pos(self, file_name=None):
        if file_name is None:
            raise Exception("OifCell: No file_name provided for vtk output. Quitting")
        n_points = len(self.mesh.points)
        n_triangles = len(self.mesh.triangles)
        output_file = open(file_name, "w")
        output_file.write("# vtk DataFile Version 3.0\n")
        output_file.write("Data\n")
        output_file.write("ASCII\n")
        output_file.write("DATASET POLYDATA\n")
        output_file.write("POINTS " + str(n_points) + " float\n")
        for p in self.mesh.points:
            coords = p.get_pos()
            output_file.write(custom_str(coords[0]) + " " + custom_str(coords[1]) + " " + custom_str(coords[2]) + "\n")
        output_file.write("TRIANGLE_STRIPS " + str(n_triangles) + " " + str(4*n_triangles) + "\n")
        for t in self.mesh.triangles:
            output_file.write("3 " + str(t.A.id) + " " + str(t.B.id) + " " + str(t.C.id) + "\n")
        output_file.close()

    def output_vtk_pos_folded(self, file_name=None):
        if file_name is None:
            raise Exception("OifCell: No file_name provided for vtk output. Quitting.")
        n_points = len(self.mesh.points)
        n_triangles = len(self.mesh.triangles)
        
        # get coordinates of the origin
        center = np.array([0.0, 0.0, 0.0])
        for p in self.mesh.points:
            center += p.get_pos()
        center /= len(self.mesh.points)
        center_folded = np.floor(center/self.cell_type.system.box_l)
        # this gives how many times the origin is folded in all three directions
        
        output_file = open(file_name, "w")
        output_file.write("# vtk DataFile Version 3.0\n")
        output_file.write("Data\n")
        output_file.write("ASCII\n")
        output_file.write("DATASET POLYDATA\n")
        output_file.write("POINTS " + str(n_points) + " float\n")
        for p in self.mesh.points:
            coords = p.get_pos() - center_folded * self.cell_type.system.box_l
            output_file.write(custom_str(coords[0]) + " " + custom_str(coords[1]) + " " + custom_str(coords[2]) + "\n")
        output_file.write("TRIANGLE_STRIPS " + str(n_triangles) + " " + str(4 * n_triangles) + "\n")
        for t in self.mesh.triangles:
            output_file.write("3 " + str(t.A.id) + " " + str(t.B.id) + " " + str(t.C.id) + "\n")
        output_file.close()

    def append_point_data_to_vtk(self, file_name=None, data_name=None, data=None, first_append=None):
        if file_name is None:
            raise Exception("OifCell: append_point_data_to_vtk: No file_name provided. Quitting.")
        if data is None:
            raise Exception("OifCell: append_point_data_to_vtk: No data provided. Quitting.")
            return
        if data_name is None:
            raise Exception("OifCell: append_point_data_to_vtk: No data_name provided. Quitting.")
        if first_append is None:
            raise Exception("OifCell: append_point_data_to_vtk: Need to know whether this is the first data list to be "
                  "appended for this file. Quitting.")
        n_points = self.get_n_nodes()
        if (len(data) != n_points):
            raise Exception("OifCell: append_point_data_to_vtk: Number of data points does not match number of mesh points. Quitting.")
        output_file = open(file_name, "a")
        if first_append is True:
            output_file.write("POINT_DATA " + str(n_points) + "\n")
        output_file.write("SCALARS " + data_name + " float 1\n")
        output_file.write("LOOKUP_TABLE default\n")
        for p in self.mesh.points:
            output_file.write(str(data[p.id]) + "\n")
        output_file.close()

    def output_raw_data(self, file_name=None, data=None):
        if file_name is None:
            raise Exception("OifCell: output_raw_data: No file_name provided. Quitting.")
        if data is None:
            raise Exception("OifCell: output_raw_data: No data provided. Quitting.")
        n_points = self.get_n_nodes()
        if (len(data) != n_points):
            raise Exception("OifCell: output_raw_data: Number of data points does not match number of mesh points. Quitting.")
        output_file = open(file_name, "w")
        for p in self.mesh.points:
            output_file.write(" ".join(map(str,data[p.id])) + "\n")
        output_file.close()

    def output_mesh_points(self, file_name=None):
        if file_name is None:
            raise Exception("OifCell: No file_name provided for mesh nodes output. Quitting.")
        output_file = open(file_name, "w")
        center = self.get_origin()
        for p in self.mesh.points:
            coords = p.get_pos() - center
            output_file.write(custom_str(coords[0]) + " " + custom_str(coords[1]) + " " + custom_str(coords[2]) + "\n")
        output_file.close()

    def set_mesh_points(self, file_name=None):
        if file_name is None:
            raise Exception("OifCell: No file_name provided for set_mesh_points. Quitting.")
        center = self.get_origin()
        n_points = self.get_n_nodes()

        in_file = open(file_name, "r")
        nodes_coord = in_file.read().split("\n")
        in_file.close()
        # removes a blank line at the end of the file if there is any:
        nodes_coord = filter(None, nodes_coord)  # here we have list of lines with triplets of strings
        if len(nodes_coord) != n_points:
            raise Exception("OifCell: Mesh nodes not set to new positions: "
                  "number of lines in the file does not equal number of Cell nodes. Quitting.")
        else:
            i = 0
            for line in nodes_coord:  # extracts coordinates from the string line
                line = line.split()
                new_position = np.array(line).astype(np.float) + center
                self.mesh.points[i].set_pos(new_position)
                i += 1

    def print_info(self):
        print("\nThe following OifCell was created: ")
        print("\t particle_mass: " + custom_str(self.particle_mass))
        print("\t particle_type: " + str(self.particle_type))
        print("\t rotate: " + str(self.rotate))
        print("\t origin: " + str(self.origin[0]) + " " + str(self.origin[1]) + " " + str(self.origin[2]))

    def elastic_forces(self, el_forces=(0, 0, 0, 0, 0, 0), f_metric=(0, 0, 0, 0, 0, 0), vtk_file=None,
                       raw_data_file=None):
        # the order of parameters in elastic_forces and in f_metric is as follows (ks, kb, kal, kag, kv, total)
        # vtk_file means that a vtk file for visualisation of elastic forces will be written
        # raw_data_file means that just the elastic forces will be written into the output file

        stretching_forces_list = []
        bending_forces_list = []
        local_area_forces_list = []
        global_area_forces_list = []
        volume_forces_list = []
        elastic_forces_list = []
        stretching_forces_norms_list = []
        bending_forces_norms_list = []
        local_area_forces_norms_list = []
        global_area_forces_norms_list = []
        volume_forces_norms_list = []
        elastic_forces_norms_list = []
        ks_f_metric = 0.0
        kb_f_metric = 0.0
        kal_f_metric = 0.0
        kag_f_metric = 0.0
        kv_f_metric = 0.0
        total_f_metric = 0.0

        for i in range(0,6):
            if (el_forces[i] != 0) and (el_forces[i] != 1):
                raise Exception("OifCell: elastic_forces: Incorrect argument. el_forces has to be a sixtuple of 0s and 1s, "
                      "specifying which elastic forces will be calculated. The order in the sixtuple is (ks, kb, "
                      "kal, kag, kv, total).")
        for i in range(0,6):
            if (f_metric[i] != 0) and (f_metric[i] != 1):
                raise Exception("OifCell: elastic_forces: Incorrect argument. f_metric has to be a sixtuple of 0s and 1s, "
                      "specifying which f_metric will be calculated. The order in the sixtuple is (ks, kb, kal, "
                      "kag, kv, total)")
        # calculation of stretching forces and f_metric
        if (el_forces[0] == 1) or (el_forces[5] == 1) or (f_metric[0] == 1) or (f_metric[5] == 1):
            # initialize list
            stretching_forces_list = []
            for p in self.mesh.points:
                stretching_forces_list.append([0.0, 0.0, 0.0])
            # calculation uses edges, but results are stored for nodes
            for e in self.mesh.edges:
                a_current_pos = e.A.get_pos()
                b_current_pos = e.B.get_pos()
                a_orig_pos = self.cell_type.mesh.points[e.A.id].get_pos()
                b_orig_pos = self.cell_type.mesh.points[e.B.id].get_pos()
                current_dist = e.length()
                orig_dist = vec_distance(a_orig_pos, b_orig_pos)
                tmp_stretching_force = oif_calc_stretching_force(self.cell_type.ks, a_current_pos, b_current_pos,
                                                             orig_dist, current_dist)
                stretching_forces_list[e.A.id] += tmp_stretching_force
                stretching_forces_list[e.B.id] -= tmp_stretching_force
            # calculation of stretching f_metric, if needed
            if f_metric[0] == 1:
                ks_f_metric = 0.0
                for p in self.mesh.points:
                    ks_f_metric += norm(stretching_forces_list[p.id])

        # calculation of bending forces and f_metric
        if (el_forces[1] == 1) or (el_forces[5] == 1) or (f_metric[1] == 1) or (f_metric[5] == 1):
            # initialize list
            bending_forces_list = []
            for p in self.mesh.points:
                bending_forces_list.append([0.0, 0.0, 0.0])
            # calculation uses bending incidences, but results are stored for nodes
            for angle in self.mesh.angles:
                a_current_pos = angle.A.get_pos()
                b_current_pos = angle.B.get_pos()
                c_current_pos = angle.C.get_pos()
                d_current_pos = angle.D.get_pos()
                a_orig_pos = self.cell_type.mesh.points[angle.A.id].get_pos()
                b_orig_pos = self.cell_type.mesh.points[angle.B.id].get_pos()
                c_orig_pos = self.cell_type.mesh.points[angle.C.id].get_pos()
                d_orig_pos = self.cell_type.mesh.points[angle.D.id].get_pos()
                current_angle = angle.size()
                orig_angle = angle_btw_triangles(a_orig_pos, b_orig_pos, c_orig_pos, d_orig_pos)
                tmp_bending_forces = oif_calc_bending_force(self.cell_type.kb, a_current_pos, b_current_pos, c_current_pos,
                                                        d_current_pos, orig_angle, current_angle)
                tmp_bending_force1 = np.array([tmp_bending_forces[0], tmp_bending_forces[1], tmp_bending_forces[2]])
                tmp_bending_force2 = np.array([tmp_bending_forces[3], tmp_bending_forces[4], tmp_bending_forces[5]])
                bending_forces_list[angle.A.id] += tmp_bending_force1
                bending_forces_list[angle.B.id] -= 0.5*tmp_bending_force1 + 0.5*tmp_bending_force2
                bending_forces_list[angle.C.id] -= 0.5*tmp_bending_force1 + 0.5*tmp_bending_force2
                bending_forces_list[angle.D.id] += tmp_bending_force2
            # calculation of bending f_metric, if needed
            if f_metric[1] == 1:
                kb_f_metric = 0.0
                for p in self.mesh.points:
                    kb_f_metric += norm(bending_forces_list[p.id])

        # calculation of local area forces and f_metric
        if (el_forces[2] == 1) or (el_forces[5] == 1) or (f_metric[2] == 1) or (f_metric[5] == 1):
            # initialize list
            local_area_forces_list = []
            for p in self.mesh.points:
                local_area_forces_list.append([0.0, 0.0, 0.0])
            # calculation uses triangles, but results are stored for nodes
            for t in self.mesh.triangles:
                a_current_pos = t.A.get_pos()
                b_current_pos = t.B.get_pos()
                c_current_pos = t.C.get_pos()
                a_orig_pos = self.cell_type.mesh.points[t.A.id].get_pos()
                b_orig_pos = self.cell_type.mesh.points[t.B.id].get_pos()
                c_orig_pos = self.cell_type.mesh.points[t.C.id].get_pos()
                current_area = t.area()
                orig_area = area_triangle(a_orig_pos, b_orig_pos, c_orig_pos)
                tmp_local_area_forces = oif_calc_local_area_force(self.cell_type.kal, a_current_pos, b_current_pos,
                                                              c_current_pos, orig_area, current_area)
                local_area_forces_list[t.A.id] += np.array([tmp_local_area_forces[0], tmp_local_area_forces[1],
                                                            tmp_local_area_forces[2]])
                local_area_forces_list[t.B.id] += np.array([tmp_local_area_forces[3], tmp_local_area_forces[4],
                                                            tmp_local_area_forces[5]])
                local_area_forces_list[t.C.id] += np.array([tmp_local_area_forces[6], tmp_local_area_forces[7],
                                                            tmp_local_area_forces[8]])

            # calculation of local area f_metric, if needed
            if f_metric[2] == 1:
                kal_f_metric = 0.0
                for p in self.mesh.points:
                    kal_f_metric += norm(local_area_forces_list[p.id])

        # calculation of global area forces and f_metric
        if (el_forces[3] == 1) or (el_forces[5] == 1) or (f_metric[3] == 1) or (f_metric[5] == 1):
            # initialize list
            global_area_forces_list = []
            for p in self.mesh.points:
                global_area_forces_list.append([0.0, 0.0, 0.0])
            # calculation uses triangles, but results are stored for nodes
            for t in self.mesh.triangles:
                a_current_pos = t.A.get_pos()
                b_current_pos = t.B.get_pos()
                c_current_pos = t.C.get_pos()
                current_surface = self.mesh.surface()
                orig_surface = self.cell_type.mesh.surface()
                tmp_global_area_forces = oif_calc_global_area_force(self.cell_type.kag, a_current_pos, b_current_pos,
                                                                c_current_pos, orig_surface, current_surface)
                global_area_forces_list[t.A.id] += np.array([tmp_global_area_forces[0], tmp_global_area_forces[1],
                                                             tmp_global_area_forces[2]])
                global_area_forces_list[t.B.id] += np.array([tmp_global_area_forces[3], tmp_global_area_forces[4],
                                                             tmp_global_area_forces[5]])
                global_area_forces_list[t.C.id] += np.array([tmp_global_area_forces[6], tmp_global_area_forces[7],
                                                             tmp_global_area_forces[8]])
            # calculation of global area f_metric, if needed
            if f_metric[3] == 1:
                kag_f_metric = 0.0
                for p in self.mesh.points:
                    kag_f_metric += norm(global_area_forces_list[p.id])

        # calculation of volume forces and f_metric
        if (el_forces[4] == 1) or (el_forces[5] == 1) or (f_metric[4] == 1) or (f_metric[5] == 1):
            # initialize list
            volume_forces_list = []
            for p in self.mesh.points:
                volume_forces_list.append([0.0, 0.0, 0.0])
            # calculation uses triangles, but results are stored for nodes
            for t in self.mesh.triangles:
                a_current_pos = t.A.get_pos()
                b_current_pos = t.B.get_pos()
                c_current_pos = t.C.get_pos()
                current_volume = self.mesh.volume()
                orig_volume = self.cell_type.mesh.volume()
                tmp_volume_force = oif_calc_volume_force(self.cell_type.kv, a_current_pos, b_current_pos, c_current_pos,
                                                     orig_volume, current_volume)
                volume_forces_list[t.A.id] += tmp_volume_force
                volume_forces_list[t.B.id] += tmp_volume_force
                volume_forces_list[t.C.id] += tmp_volume_force
            # calculation of volume f_metric, if needed
            if f_metric[4] == 1:
                kv_f_metric = 0.0
                for p in self.mesh.points:
                    kv_f_metric += norm(volume_forces_list[p.id])

        # calculation of total elastic forces and f_metric
        if (el_forces[5] == 1) or (f_metric[5] == 1):
            elastic_forces_list = []
            for p in self.mesh.points:
                total_elastic_forces = stretching_forces_list[p.id] + bending_forces_list[p.id] + \
                                       local_area_forces_list[p.id] + global_area_forces_list[p.id] + \
                                       volume_forces_list[p.id]
                elastic_forces_list.append(total_elastic_forces)
            # calculation of total f_metric, if needed
            if f_metric[5] == 1:
                total_f_metric = 0.0
                for p in self.mesh.points:
                    total_f_metric += norm(elastic_forces_list[p.id])

        # calculate norms of resulting forces
        if (el_forces[0] + el_forces[1] + el_forces[2] + el_forces[3] + el_forces[4] + el_forces[5]) != 0:
            if el_forces[0] == 1:
                stretching_forces_norms_list = []
                for p in self.mesh.points:
                    stretching_forces_norms_list.append(norm(stretching_forces_list[p.id]))
            if el_forces[1] == 1:
                bending_forces_norms_list = []
                for p in self.mesh.points:
                    bending_forces_norms_list.append(norm(bending_forces_list[p.id]))
            if el_forces[2] == 1:
                local_area_forces_norms_list = []
                for p in self.mesh.points:
                    local_area_forces_norms_list.append(norm(local_area_forces_list[p.id]))
            if el_forces[3] == 1:
                global_area_forces_norms_list = []
                for p in self.mesh.points:
                    global_area_forces_norms_list.append(norm(global_area_forces_list[p.id]))
            if el_forces[4] == 1:
                volume_forces_norms_list = []
                for p in self.mesh.points:
                    volume_forces_norms_list.append(norm(volume_forces_list[p.id]))
            if el_forces[5] == 1:
                elastic_forces_norms_list = []
                for p in self.mesh.points:
                    elastic_forces_norms_list.append(norm(elastic_forces_list[p.id]))

        # output vtk (folded)
        if vtk_file is not None:
            if el_forces == (0, 0, 0, 0, 0, 0):
                raise Exception("OifCell: elastic_forces: The option elastic_forces was not used. "
                      "Nothing to output to vtk file.")
            self.output_vtk_pos_folded(vtk_file)
            first = True
            if el_forces[0] == 1:
                self.append_point_data_to_vtk(file_name=vtk_file, data_name="ks_f_metric",
                                              data=stretching_forces_norms_list, first_append=first)
                first = False
            if el_forces[1] == 1:
                self.append_point_data_to_vtk(file_name=vtk_file, data_name="kb_f_metric",
                                              data=bending_forces_norms_list, first_append=first)
                first = False
            if el_forces[2] == 1:
                self.append_point_data_to_vtk(file_name=vtk_file, data_name="kal_f_metric",
                                              data=local_area_forces_norms_list, first_append=first)
                first = False
            if el_forces[3] == 1:
                self.append_point_data_to_vtk(file_name=vtk_file, data_name="kag_f_metric",
                                              data=global_area_forces_norms_list, first_append=first)
                first = False
            if el_forces[4] == 1:
                self.append_point_data_to_vtk(file_name=vtk_file, data_name="kav_f_metric",
                                              data=volume_forces_norms_list, first_append=first)
                first = False
            if el_forces[5] == 1:
                self.append_point_data_to_vtk(file_name=vtk_file, data_name="total_f_metric",
                                              data=elastic_forces_norms_list, first_append=first)
                first = False

        # output raw data
        if raw_data_file is not None:
            if (el_forces[0] + el_forces[1] + el_forces[2] + el_forces[3] + el_forces[4] + el_forces[5]) != 1:
                raise Exception("OifCell: elastic_forces: Only one type of elastic forces can be written into one "
                      "raw_data_file. If you need several, please call OifCell.elastic_forces multiple times - "
                      "once per elastic force.")
            if el_forces[0] == 1:
                self.output_raw_data(file_name=raw_data_file, data=stretching_forces_list)
            if el_forces[1] == 1:
                self.output_raw_data(file_name=raw_data_file, data=bending_forces_list)
            if el_forces[2] == 1:
                self.output_raw_data(file_name=raw_data_file, data=local_area_forces_list)
            if el_forces[3] == 1:
                self.output_raw_data(file_name=raw_data_file, data=global_area_forces_list)
            if el_forces[4] == 1:
                self.output_raw_data(file_name=raw_data_file, data=volume_forces_list)
            if el_forces[5] == 1:
                self.output_raw_data(file_name=raw_data_file, data=elastic_forces_list)

        # return f_metric
        if f_metric[0] + f_metric[1] + f_metric[2] + f_metric[3] + f_metric[4] + f_metric[5] > 0:
            results = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            if f_metric[0] == 1:
                results[0] = ks_f_metric
            if f_metric[1] == 1:
                results[1] = kb_f_metric
            if f_metric[2] == 1:
                results[2] = kal_f_metric
            if f_metric[3] == 1:
                results[3] = kag_f_metric
            if f_metric[4] == 1:
                results[4] = kv_f_metric
            if f_metric[5] == 1:
                results[5] = total_f_metric
            return results
        else:
            return 0
