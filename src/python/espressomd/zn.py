#
# Copyright (C) 2024-2026 The ESPResSo project
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

import subprocess
import numpy as np
import zndraw
import zndraw.utils
import zndraw.materials
import zndraw.geometries
import espressomd
import secrets
import time
import urllib.parse
import requests
import typing
import scipy.spatial.transform

from espressomd.plugins import ase


class LBField:
    """
    Convert the ESPResSo lattice-Boltzmann field to a vector field for visualization. Samples
    the field at a given step size and offset over the lattice nodes.

    Parameters
    ----------
    system : :class:`~espressomd.system.System`
        ESPResSo system
    step_x : :obj:`int`, optional
        Step size in x direction, by default 1
    step_y : :obj:`int`, optional
        Step size in y direction, by default 1
    step_z : :obj:`int`, optional
        Step size in z direction, by default 1
    offset_x : :obj:`int`, optional
        Offset in x direction, by default 0
    offset_y : :obj:`int`, optional
        Offset in y direction, by default 0
    offset_z : :obj:`int`, optional
        Offset in z direction, by default 0
    scale : :obj:`float`, optional
        Scale the velocity vectors, by default 1.0
    arrow_config : :obj:`dict`, optional
        Configuration for the arrows, by default None and then uses the default configuration:

        'colormap': [[-0.5, 0.9, 0.5], [-0.0, 0.9, 0.5]]
            HSL colormap for the arrows, where the first value is the minimum value and the second value is the maximum value.
        'normalize': True
            Normalize the colormap to the maximum value each frame
        'colorrange': [0, 1]
            Range of the colormap, only used if normalize is False
        'scale_vector_thickness': True
            Scale the thickness of the arrows with the velocity
        'opacity': 1.0
            Opacity of the arrows
    """

    def __init__(self, system: espressomd.system.System,
                 step_x: int = 1,
                 step_y: int = 1,
                 step_z: int = 1,
                 offset_x: int = 0,
                 offset_y: int = 0,
                 offset_z: int = 0,
                 scale: float = 1.0,
                 arrow_config: dict = None):
        self.step_x = step_x
        self.step_y = step_y
        self.step_z = step_z
        self.offset_x = offset_x
        self.offset_y = offset_y
        self.offset_z = offset_z
        self.scale = scale
        self.box = system.box_l

        if system.lb is None:
            raise ValueError("System does not have a lattice-Boltzmann solver")
        self.lbf = system.lb
        self.agrid = system.lb.agrid
        self.arrow_config = arrow_config

    def _build_grid(self):
        x = np.arange(self.offset_x + self.agrid / 2,
                      self.box[0], self.step_x * self.agrid)
        y = np.arange(self.offset_y + self.agrid / 2,
                      self.box[1], self.step_y * self.agrid)
        z = np.arange(self.offset_z + self.agrid / 2,
                      self.box[2], self.step_z * self.agrid)

        origins = np.array(np.meshgrid(x, y, z)).T.reshape(-1, 3)

        return origins

    def _get_velocities(self):
        velocities = self.lbf[:, :, :].velocity
        velocities = velocities[::self.step_x, ::self.step_y, ::self.step_z]
        velocities = self.scale * velocities
        velocities = np.swapaxes(velocities, 0, 3)
        velocities = np.swapaxes(velocities, 2, 3)
        velocities = velocities.T.reshape(-1, 3)

        return velocities

    def __call__(self):
        origins = self._build_grid()
        velocities = self._get_velocities()
        vector_field = np.stack([origins, velocities], axis=1)

        return vector_field


class VectorField:
    """
    Give an array of origins and vectors to create a vector field for visualization.
    The vectorfield is updated every time it is called. Both origins and vectors must have the same shape.
    Both origins and vectors must be 3D numpy arrays in the shape of ``(n, m, 3)``,
    where ``n`` is the number of frames the field has and ``m`` is the number of vectors.
    The number of frames ``n`` must larger or equal to the number of times the update function will be called in the update loop.

    Parameters
    ----------
    origins : (n, m, 3) array_like of :obj:`float`
        Array of origins for the vectors
    vectors : (n, m, 3) array_like of :obj:`float`
        Array of vectors
    scale : :obj:`float`, optional
        Scale the vectors, by default 1
    arrow_config : :obj:`dict`, optional
        Configuration for the arrows, by default ``None`` and then uses the default configuration:

        'colormap': [[0.5, 0.9, 0.5], [0.0, 0.9, 0.5]]
            HSL colormap for the arrows, where the first value is the minimum value and the second value is the maximum value.
        'normalize': True
            Normalize the colormap to the maximum value each frame
        'opacity': 1.0
            Opacity of the arrows
    """

    def __init__(self, origins: np.ndarray, vectors: np.ndarray,
                 scale: float = 1, arrow_config: dict = None):
        self.origins = origins
        self.vectors = vectors
        self.scale = scale
        self.frame_count = 0
        self.arrow_config = arrow_config

    def __call__(self):
        if self.origins.shape != self.vectors.shape:
            raise ValueError("Origins and vectors must have the same shape")

        origins = self.origins[self.frame_count]
        vectors = origins + self.vectors[self.frame_count]
        self.frame_count += 1
        origins = origins.reshape(-1, 3)
        vectors = vectors.reshape(-1, 3)
        vectors = self.scale * vectors
        vector_field = np.stack([origins, vectors], axis=1)

        return vector_field


class Visualizer():
    """
    Visualizer for ESPResSo simulations using ZnDraw.

    Main component of the visualizer is the ZnDraw server, which is started as a subprocess.
    The ZnDraw client is used to communicate with the server and send the visualized data.
    The visualizer uploads a new frame to the server every time the update method is called.

    Parameters
    ----------
    system : :class:`~espressomd.system.System`
        ESPResSo system to visualize
    port : :obj:`int`, optional
        Port for the ZnDraw server, by default 1234, if taken, the next available port is used
    token : :obj:`str`, optional
        Token for the ZnDraw server, by default a random token is generated
    folded : :obj:`bool`, optional
        Fold the positions of the particles into the simulation box, by default True
    colors : :obj:`dict`, optional
        Dictionary containing color type mapping for the particles, by default all particles are white
    radii : :obj:`dict`, optional
        Dictionary containing radii type mapping for the particles, by default all particles have a radius of 0.5
    bonds : :obj:`bool`, optional
        Draw bonds between particles, by default ``False``
    jupyter : :obj:`bool`, optional
        Show the visualizer in a Jupyter notebook, by default True
    vector_field : :class:`~espressomd.zn.VectorField` or :class:`~espressomd.zn.LBField`, optional
        Vector field to visualize, by default ``None``

    """

    SERVER_PORT = None

    def __init__(self,
                 system: espressomd.system.System,
                 port: int = 1234,
                 token: str = None,
                 folded: bool = True,
                 colors: dict = None,
                 radii: dict = None,
                 bonds: bool = False,
                 forces: bool = False,
                 jupyter: bool = True,
                 vector_field: typing.Union[VectorField, LBField] = None,
                 ):

        self.system = system
        self.params = {
            "folded": folded,
            "colors": colors,
            "radii": radii,
            "bonds": bonds,
            "forces": forces,
            "vector_field": vector_field,
        }

        self.url = "http://127.0.0.1"
        self.frame_count = 0
        if token is None:
            self.token = secrets.token_hex(4)
        else:
            self.token = token

        # A server is started in a subprocess, and we have to wait for it
        if self.SERVER_PORT is None:
            print("Starting ZnDraw server, this may take a few seconds")
            Visualizer.SERVER_PORT = port
            self.server = subprocess.Popen(
                ["zndraw", "--no-browser", f"--port={self.SERVER_PORT}"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL)

        # Check that the server is up
        request_deadline = time.monotonic() + 30
        while time.monotonic() < request_deadline:
            try:
                r = requests.get(
                    f"{self.url}:{self.SERVER_PORT}/health", timeout=1)
                if r.status_code == 200:
                    break
            except requests.RequestException:
                pass
            time.sleep(0.5)
        else:
            raise RuntimeError(
                "ZnDraw server did not start within the expected time")

        self._start_zndraw()

        if vector_field is not None:
            self.arrow_config = {'colormap': [[0.5, 0.9, 0.5], [0.0, 0.9, 0.5]],
                                 'normalize': True,
                                 'opacity': 1.0}

            if vector_field.arrow_config is not None:
                for key, value in vector_field.arrow_config.items():
                    if key not in self.arrow_config:
                        raise ValueError(f"Invalid key {key} in arrow_config")
                    self.arrow_config[key] = value

        if jupyter:
            self._show_jupyter()
        else:
            raise NotImplementedError(
                "Only Jupyter notebook is supported at the moment")

        # Wait until the session is ready
        request_deadline = time.monotonic() + 30
        while time.monotonic() < request_deadline:
            if len(self.zndraw.sessions) > 0:
                break
            time.sleep(0.5)
        else:
            raise RuntimeError(
                "ZnDraw session did not start within the expected time")

    def _start_zndraw(self):
        """
        Start the ZnDraw client and connect to the server
        """
        url = f"{self.url}:{self.SERVER_PORT}"
        self.zndraw = zndraw.ZnDraw(url=url, room=self.token)
        parsed_url = urllib.parse.urlparse(
            f"{self.zndraw.url}/rooms/{self.zndraw.room}")
        self.address = parsed_url._replace(scheme="http").geturl()

    def _show_jupyter(self):
        """
        Show the visualizer in a Jupyter notebook
        """
        from IPython.display import IFrame, display
        print(f"Showing ZnDraw at {self.address}")
        display(IFrame(src=self.address, width="100%", height="700px"))

    def update(self):
        """
        Update the visualizer with the current state of the system
        """
        particles = self.system.part.all()
        all_types = particles.type
        self.system.visualizer_params = self.params
        ase_interf = ase.ASEInterface(
            system=self.system,
            particle_slice=self.system.part.all(),
            use_folded_positions=self.params["folded"],
            export_forces=self.params["forces"])
        data = ase_interf.atoms
        if self.params["colors"] is not None:
            data.arrays["colors"] = [self.params["colors"].get(
                z, "white") for z in all_types]
        if self.params["radii"] is not None:
            data.arrays["radii"] = [self.params["radii"].get(
                z, 0.5) for z in all_types]

        if self.params["bonds"]:
            id_mapping = {p.id: idx for idx, p in enumerate(particles)}
            bonds = []
            for p in particles:
                if not p.bonds:
                    continue
                for bond in p.bonds:
                    if len(bond) == 4:
                        bonds.extend([
                            [id_mapping[p.id], id_mapping[bond[1]], 1],
                            [id_mapping[p.id], id_mapping[bond[2]], 1],
                            [id_mapping[bond[2]], id_mapping[bond[3]], 1]
                        ])
                    else:
                        bonds.extend(
                            [[id_mapping[p.id], id_mapping[bond_partner], 1] for bond_partner in bond[1:]])

            if self.params["folded"]:
                # add ghost particles and bonds for PBC crossing bonds
                bonds = self._handle_pbc_bonds(bonds=bonds, ase_data=data)
            data.info["connectivity"] = bonds

        if self.frame_count == 0:
            x, y, z = self.system.box_l / 2
            z_dist = max([1.5 * y, 1.5 * x, 1.5 * z])

            self.zndraw.geometries["camera"] = zndraw.geometries.Camera(
                position=(x, y, z_dist),
                target=(x, y, z),
            )

            # activate the camera
            for key in self.zndraw.sessions.keys():
                self.zndraw.api.set_active_camera(key, "camera")

            if self.params["forces"]:
                self.zndraw.geometries["force_arrows"] = zndraw.geometries.Arrow(
                    position="arrays.positions",
                    direction="arrays.forces",
                    color=["gray"] if self.params["colors"] else "arrays.colors",
                    radius=4 *
                    max(self.params["radii"].values()
                        ) if self.params["radii"] else 1.0,
                    opacity=0.8,
                    hovering=zndraw.geometries.InteractionSettings(
                        enabled=False),
                    selecting=zndraw.geometries.InteractionSettings(
                        enabled=False),
                )

        if self.params["vector_field"] is not None:
            field_data = self.params["vector_field"]()
            data.info["vector_positions"] = field_data[:, 0]
            vector_directions = np.copy(field_data[:, 1])
            magnitudes = np.linalg.norm(
                vector_directions, axis=1, keepdims=True)
            if self.arrow_config['normalize']:
                mask = magnitudes.flatten() > 0.
                vector_directions[mask] /= magnitudes[mask]
            data.info["vector_directions"] = vector_directions

            colormap = np.array(self.arrow_config['colormap'])
            max_magnitude = np.max(magnitudes)
            min_magnitude = np.min(magnitudes)
            vector_color = colormap[0] + (magnitudes - min_magnitude) / (
                max_magnitude - min_magnitude) * (colormap[1] - colormap[0])
            vector_color = _convert_array_to_color(vector_color)
            data.info["vector_color"] = vector_color

            for key, value in self.arrow_config.items():
                data.info[key] = value
            if "vector_field" not in self.zndraw.geometries.keys():
                shape_options = {"hovering": zndraw.geometries.InteractionSettings(enabled=False),
                                 "selecting": zndraw.geometries.InteractionSettings(enabled=False)}
                self.zndraw.geometries["vector_field"] = zndraw.geometries.Arrow(
                    position="info.vector_positions",
                    direction="info.vector_directions",
                    color="info.vector_color",
                    **shape_options)

        # Catch when the server is initializing an empty frame
        # len(self.zndraw) is a expensive socket call, so we try to avoid it
        if self.frame_count != 0 or len(self.zndraw) == 0:
            self.zndraw.append(data)
        else:
            self.zndraw[0] = data

        self.frame_count += 1

    def register_setting(self, cls, **kwargs):
        self.zndraw.register_modifier(cls, **kwargs)

    def draw_constraints(self, shapes: list, **shape_options):
        """
        Draw constraints on the visualizer
        """
        if not isinstance(shapes, list):
            raise ValueError("Constraints must be given in a list")

        if "hovering" not in shape_options:
            shape_options["hovering"] = zndraw.geometries.InteractionSettings(
                enabled=False)
        if "selecting" not in shape_options:
            shape_options["selecting"] = zndraw.geometries.InteractionSettings(
                enabled=False)

        for shape in shapes:
            shape_type = shape.__class__.__name__

            if shape_type == "Sphere":
                center = tuple(shape.center)
                radius = shape.radius
                key = f"{shape_type}_{center}_{radius}"

                self.zndraw.geometries[key] = zndraw.geometries.Sphere(
                    position=[tuple(center)], radius=[radius], **shape_options)

            elif shape_type == "Wall":
                dist = shape.dist
                normal = np.array(shape.normal)

                position = dist * normal
                helper = WallIntersection(
                    plane_normal=normal, plane_point=position, box_l=self.system.box_l)
                corners = helper.get_intersections()

                base_position = np.copy(corners[0])
                corners -= base_position

                vertices, euler_angles = _corners_to_shape_geometry(corners)

                key = f"{shape_type}_{dist}_{normal}"
                self.zndraw.geometries[key] = zndraw.geometries.Shape(
                    position=[tuple(base_position)],
                    rotation=[tuple(euler_angles)],
                    vertices=vertices, **shape_options)

            elif shape_type == "Rhomboid":
                vecs = np.array([shape.a, shape.b, shape.c])
                corner_base = shape.corner
                for direction in range(3):
                    for side in range(2):
                        vec1 = vecs[(direction + 1) % 3]
                        vec2 = vecs[(direction + 2) % 3]
                        corner = np.copy(corner_base)
                        if (side == 1):
                            corner += vecs[direction]
                        corners = np.array([corner,
                                            corner + vec1,
                                            corner + vec1 + vec2,
                                            corner + vec2])
                        base_position = np.copy(corners[0])
                        corners -= base_position
                        vertices, euler_angles = _corners_to_shape_geometry(
                            corners, False)
                        key = f"{shape_type}_{corner_base}_{vecs}_{direction * 2 + side}"  # nopep8
                        self.zndraw.geometries[key] = zndraw.geometries.Shape(
                            position=[tuple(base_position)],
                            rotation=[tuple(euler_angles)],
                            vertices=vertices, **shape_options)

            else:
                raise NotImplementedError(
                    f"Shape of type {shape_type} isn't available in ZnDraw")

    def _handle_pbc_bonds(self, bonds, ase_data):
        box_l = self.system.box_l
        positions = ase_data.arrays['positions']
        max_index = len(positions)
        colors = ase_data.arrays['colors']
        min_radii = np.min(ase_data.arrays['radii'])
        numbers = ase_data.arrays['numbers']
        bonds_to_add = []
        ghost_positions = []
        ghost_colors = []
        ghost_numbers = []

        for index, b in reversed(list(enumerate(bonds))):
            pos_a = positions[b[0]]
            pos_b = positions[b[1]]

            distance_vec = pos_b - pos_a

            if np.all(np.abs(distance_vec) < 0.5 * box_l):
                continue

            # bond is crossing the PBC
            bonds.pop(index)

            # get the minimum-image bond vector
            distance_vec -= np.rint(distance_vec / box_l) * box_l

            # find the wall intersection point
            best_distance = np.inf
            for i in range(3):
                if distance_vec[i] == 0:
                    continue  # corner case: bond is parallel to a face
                elif distance_vec[i] > 0:
                    p0_i = box_l[i]
                else:
                    p0_i = 0
                # calculate the length of the bond that is inside the box using
                # an optimized version of `np.dot(p0 - x0, n) / np.dot(dx, n)`
                # where the dot products decay to an array access, because `n`
                # is always a unit vector in rectangular boxes and its sign
                # is canceled out by the division
                d = (p0_i - pos_a[i]) / distance_vec[i]
                if d < best_distance:
                    best_distance = d

            wall_intersection_a = pos_a + best_distance * distance_vec
            wall_intersection_b = pos_b - (1 - best_distance) * distance_vec

            ghost_positions.extend([wall_intersection_a, wall_intersection_b])
            ghost_colors.extend([colors[b[0]], colors[b[1]]])
            ghost_numbers.extend([numbers[b[0]], numbers[b[1]]])
            bonds_to_add.extend(
                [(b[0], max_index, 1), (b[1], max_index + 1, 1)])

            max_index += 2

        if not ghost_positions:
            return bonds

        # add ghost particles
        ase_data.arrays['positions'] = np.vstack([positions, ghost_positions])
        ase_data.arrays['colors'] = np.hstack(
            [ase_data.arrays['colors'], ghost_colors])
        ase_data.arrays['radii'] = np.hstack(
            [ase_data.arrays['radii'], [1e-6 * min_radii] * len(ghost_positions)])
        ase_data.arrays['numbers'] = np.hstack([numbers, ghost_numbers])
        ase_data.arrays['forces'] = np.vstack(
            [ase_data.arrays['forces'], np.zeros((len(ghost_positions), 3))])

        bonds.extend(bonds_to_add)
        return bonds


class WallIntersection:
    """
    Simple helper to calculate all Box edges that intersect with a plane.
    """

    def __init__(self, plane_point, plane_normal, box_l):
        self.plane_point = plane_point
        self.plane_normal = plane_normal / np.linalg.norm(plane_normal)
        self.box_l = box_l

        # Create 8 vertices of the bounding box
        self.vertices = np.array([
            [0, 0, 0],
            [0, 0, box_l[2]],
            [0, box_l[1], 0],
            [0, box_l[1], box_l[2]],
            [box_l[0], 0, 0],
            [box_l[0], 0, box_l[2]],
            [box_l[0], box_l[1], 0],
            [box_l[0], box_l[1], box_l[2]]
        ])

        self.edges = [
            (self.vertices[0], self.vertices[1]),
            (self.vertices[0], self.vertices[2]),
            (self.vertices[0], self.vertices[4]),
            (self.vertices[1], self.vertices[3]),
            (self.vertices[1], self.vertices[5]),
            (self.vertices[2], self.vertices[3]),
            (self.vertices[2], self.vertices[6]),
            (self.vertices[3], self.vertices[7]),
            (self.vertices[4], self.vertices[5]),
            (self.vertices[4], self.vertices[6]),
            (self.vertices[5], self.vertices[7]),
            (self.vertices[6], self.vertices[7])
        ]

    def plane_intersection_with_line(self, line_point1, line_point2):
        # Calculate the intersection point of a line and a plane
        line_dir = line_point2 - line_point1
        denom = np.dot(self.plane_normal, line_dir)

        if np.abs(denom) > 1e-6:  # Avoid division by zero
            t = np.dot(self.plane_normal,
                       (self.plane_point - line_point1)) / denom
            if 0 <= t <= 1:
                return line_point1 + t * line_dir
        return None

    def get_intersections(self):
        intersections = []

        for edge in self.edges:
            intersection = self.plane_intersection_with_line(edge[0], edge[1])
            if intersection is not None:
                intersections.append(intersection)

        return np.array(intersections)


def _corners_to_shape_geometry(corners, sort=True):
    """
    Calculate the vertices and Euler angles of a rhombus defined by its corners.
    """
    unit_z = np.array([0., 0., 1.])
    v1 = corners[1] - corners[0]
    v2 = corners[-1] - corners[0]
    normal = np.cross(v2, v1)
    normal = normal / np.linalg.norm(normal)

    rot, _ = scipy.spatial.transform.Rotation.align_vectors(
        [normal], [unit_z])

    rot_matix = np.linalg.inv(rot.as_matrix())
    vertices = np.column_stack([
        corners @ rot_matix[0, :],
        corners @ rot_matix[1, :]
    ])

    if sort:
        angles = np.arctan2(vertices[:, 1], vertices[:, 0])
        sorted_indices = np.argsort(angles)
        sorted_vertices = vertices[:][sorted_indices]
    else:
        sorted_vertices = vertices

    euler_angles = rot.as_euler('xyz')
    # invert the z-axis (different coordinate system)
    euler_angles[2] *= -1

    return sorted_vertices, euler_angles


def _convert_array_to_color(array):
    """
    Converts a (n,3) array to a list of color strings.
    """
    if array.shape[-1] != 3:
        raise NotImplementedError(
            f"The color array should be 3 dimensional")
    hex_array = array * 255
    hex_array = np.clip(hex_array, 0, 255)
    hex_array = hex_array.astype(np.uint8)
    hex_array = hex_array.reshape((-1, hex_array.shape[-1]))
    hex_strings = [f"#{r:02X}{g:02X}{b:02X}" for r, g, b in hex_array]
    return hex_strings
