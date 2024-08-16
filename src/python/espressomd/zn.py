#
# Copyright (C) 2024 The ESPResSo project
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
import zndraw.zndraw
import zndraw.utils
import zndraw.draw
import znsocket
import znjson
import espressomd
import secrets
import time
import urllib.parse
import typing as t
import scipy.spatial.transform


# Standard colors
color_dict = {"black": "#303030",
              "red": "#e6194B",
              "green": "#3cb44b",
              "yellow": "#ffe119",
              "blue": "#4363d8",
              "orange": "#f58231",
              "purple": "#911eb4",
              "cyan": "#42d4f4",
              "magenta": "#f032e6",
              "lime": "#bfef45",
              "brown": "#9A6324",
              "grey": "#a9a9a9",
              "white": "#f0f0f0"}


class EspressoConverter(znjson.ConverterBase):
    """
    Converter for ESPResSo systems to ASEDict
    """
    level = 100
    representation = "ase.Atoms"
    instance = espressomd.system.System

    def encode(self, system) -> zndraw.utils.ASEDict:
        self.system = system
        self.particles = self.system.part.all()
        self.num_particles = len(self.particles)
        self.params = system.visualizer_params

        self.numbers = self.num_particles * [1]

        if self.params["folded"] is True:
            self.positions = self.particles.pos_folded
        else:
            self.positions = self.particles.pos

        if self.params["colors"] is None:
            self.colors = self.get_default_colors()
        else:
            self.colors = self.set_colors(self.params["colors"])

        if self.params["radii"] is None:
            self.radii = self.get_default_radii()
        else:
            self.radii = self.set_radii(self.params["radii"])

        if self.params["bonds"] is True:
            bonds = self.get_bonds()
        else:
            bonds = []

        arrays = {
            "colors": self.colors,
            "radii": self.radii,
        }
        cell = [[system.box_l[0], 0, 0],
                [0, system.box_l[1], 0],
                [0, 0, system.box_l[2]]]
        pbc = system.periodicity
        calc = None
        info = {}

        if self.params["vector_field"] is not None:
            vectors = self.params["vector_field"]()
        else:
            vectors = []

        return zndraw.utils.ASEDict(
            numbers=self.numbers,
            positions=self.positions.tolist(),
            connectivity=bonds,
            arrays=arrays,
            info=info,
            calc=calc,
            pbc=pbc.tolist(),
            cell=cell,
            vectors=vectors,
        )

    def decode(self, value):
        value = None
        return value

    def get_default_colors(self):
        return [color_dict["white"]] * self.num_particles

    def get_default_radii(self):
        return [0.5] * self.num_particles

    def set_colors(self, colors):
        color_list = list()
        for p in self.particles:
            color = colors[p.type]
            # if color starts with #, assume it is a hex color
            if color.startswith("#"):
                color_list.append(color)
            else:
                if color not in color_dict:
                    raise ValueError(
                        f"Color {color} not found in color dictionary")
                color_list.append(color_dict[color])
        return color_list

    def set_radii(self, radii):
        radius_list = list()
        for p in self.particles:
            radius_list.append(radii[p.type])
        return radius_list

    def get_bonds(self):
        bonds = []
        for p in self.particles:
            if not p.bonds:
                continue
            for bond in p.bonds:
                if len(bond) == 4:
                    bonds.append([p.id, bond[1], 1])
                    bonds.append([p.id, bond[2], 1])
                    bonds.append([bond[2], bond[3], 1])
                else:
                    for bond_partner in bond[1:]:
                        bonds.append([p.id, bond_partner, 1])

        self.process_bonds(bonds)

        return bonds

    def process_bonds(self, bonds):
        half_box_l = 0.5 * self.system.box_l
        num_part = len(self.positions)
        bonds_to_remove = []
        bonds_to_add = []

        for b in bonds:
            try:
                if self.params["folded"] is True:
                    x_a = self.system.part.by_id(b[0]).pos_folded
                    x_b = self.system.part.by_id(b[1]).pos_folded
                else:
                    x_a = self.system.part.by_id(b[0]).pos
                    x_b = self.system.part.by_id(b[1]).pos
            except Exception:
                bonds_to_remove.append(b)
                continue

            dx = x_b - x_a

            if np.all(np.abs(dx) < half_box_l):
                continue

            if self.params["folded"] is False:
                bonds_to_remove.append(b)
                continue

            d = self.cut_bond(x_a, dx)
            if d is np.inf:
                bonds_to_remove.append(b)
                continue

            s_a = x_a + 0.8 * dx
            s_b = x_b - 0.8 * dx

            bonds_to_remove.append(b)

            self.add_ghost_particle(pos=s_a, color=self.colors[b[0]])
            bonds_to_add.append([b[0], num_part, 1])

            self.add_ghost_particle(pos=s_b, color=self.colors[b[1]])
            bonds_to_add.append([b[1], num_part + 1, 1])
            num_part += 2

        for b in bonds_to_remove:
            bonds.remove(b)

        bonds.extend(bonds_to_add)

    def cut_bond(self, x_a, dx):
        if np.dot(dx, dx) < 1e-9:
            return np.inf
        shift = np.rint(dx / self.system.box_l)
        dx -= shift * self.system.box_l
        best_d = np.inf
        for i in range(3):
            if dx[i] == 0:
                continue
            elif dx[i] > 0:
                p0_i = self.system.box_l[i]
            else:
                p0_i = 0

            d = (p0_i - x_a[i]) / dx[i]
            if d < best_d:
                best_d = d
        return best_d

    def add_ghost_particle(self, pos, color):
        self.positions = np.vstack([self.positions, pos])
        self.radii.append(1e-6 * min(self.radii))
        self.colors.append(color)
        self.numbers.append(2)


znjson.config.register(EspressoConverter)


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
        velocities = origins + velocities
        vector_field = np.stack([origins, velocities], axis=1)

        return vector_field.tolist()


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

        return vector_field.tolist()


class Visualizer():
    """
    Visualizer for ESPResSo simulations using ZnDraw.

    Main component of the visualizer is the ZnDraw server, which is started as a subprocess.
    The ZnDraw client is used to communicate with the server and send the visualized data.
    The visualized data is encoded using the :class:`EspressoConverter`,
    which converts the ESPResSo system to an typed dict.
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
    SOCKET_PORT = None

    def __init__(self,
                 system: espressomd.system.System = None,
                 port: int = 1234,
                 token: str = None,
                 folded: bool = True,
                 colors: dict = None,
                 radii: dict = None,
                 bonds: bool = False,
                 jupyter: bool = True,
                 vector_field: t.Union[VectorField, LBField] = None,
                 ):

        self.system = system
        self.params = {
            "folded": folded,
            "colors": colors,
            "radii": radii,
            "bonds": bonds,
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
            self.port = port
            self._start_server()
            time.sleep(10)

        self._start_zndraw()
        time.sleep(1)

        if vector_field is not None:
            self.arrow_config = {'colormap': [[-0.5, 0.9, 0.5], [-0.0, 0.9, 0.5]],
                                 'normalize': True,
                                 'colorrange': [0, 1],
                                 'scale_vector_thickness': True,
                                 'opacity': 1.0}

            if vector_field.arrow_config is not None:
                for key, value in vector_field.arrow_config.items():
                    if key not in self.arrow_config:
                        raise ValueError(f"Invalid key {key} in arrow_config")
                    self.arrow_config[key] = value

        if self.params["bonds"] and not self.params["folded"]:
            print(
                "Warning: Unfolded positions may result in incorrect bond visualization")

        if jupyter:
            self._show_jupyter()
        else:
            raise NotImplementedError(
                "Only Jupyter notebook is supported at the moment")

    def _start_server(self):
        """
        Start the ZnDraw server through a subprocess
        """
        self.socket_port = zndraw.utils.get_port(default=6374)

        Visualizer.SERVER_PORT = self.port
        Visualizer.SOCKET_PORT = self.socket_port

        self.server = subprocess.Popen(["zndraw", "--no-browser", f"--port={self.port}", f"--storage-port={self.socket_port}"],
                                       stdout=subprocess.DEVNULL,
                                       stderr=subprocess.DEVNULL
                                       )

    def _start_zndraw(self):
        """
        Start the ZnDraw client and connect to the server
        """
        config = zndraw.zndraw.TimeoutConfig(
            connection=10,
            modifier=0.25,
            between_calls=0.1,
            emit_retries=3,
            call_retries=1,
            connect_retries=3,
        )
        while True:
            try:
                self.r = znsocket.Client(
                    address=f"{self.url}:{self.SOCKET_PORT}")
                break
            except BaseException:
                time.sleep(0.5)

        url = f"{self.url}:{self.SERVER_PORT}"
        self.zndraw = zndraw.zndraw.ZnDrawLocal(
            r=self.r, url=url, token=self.token, timeout=config)
        parsed_url = urllib.parse.urlparse(
            f"{self.zndraw.url}/token/{self.zndraw.token}")
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
        self.system.visualizer_params = self.params

        data = znjson.dumps(
            self.system, cls=znjson.ZnEncoder.from_converters(
                [EspressoConverter])
        )

        # Catch when the server is initializing an empty frame
        # len(self.zndraw) is a expensive socket call, so we try to avoid it
        if self.frame_count != 0 or len(self.zndraw) == 0:
            self.zndraw.append(data)
        else:
            self.zndraw.__setitem__(0, data)

        if self.frame_count == 0:
            self.zndraw.socket.sleep(1)

            x = self.system.box_l[0] / 2
            y = self.system.box_l[1] / 2
            z = self.system.box_l[2] / 2

            z_dist = max([1.5 * y, 1.5 * x, 1.5 * z])

            self.zndraw.camera = {'position': [
                x, y, z_dist], 'target': [x, y, z]}
            self.zndraw.config.scene.frame_update = False

            if self.params["vector_field"] is not None:
                for key, value in self.arrow_config.items():
                    setattr(self.zndraw.config.arrows, key, value)

        self.frame_count += 1

    def draw_constraints(self, shapes: list):
        """
        Draw constraints on the visualizer
        """
        if not isinstance(shapes, list):
            raise ValueError("Constraints must be given in a list")

        objects = []

        for shape in shapes:

            shape_type = shape.__class__.__name__

            mat = zndraw.draw.Material(color="#b0b0b0", opacity=0.8)

            if shape_type == "Cylinder":
                center = shape.center
                axis = shape.axis
                length = shape.length
                radius = shape.radius

                rotation_angles = zndraw.utils.direction_to_euler(
                    axis, roll=np.pi / 2)

                objects.append(zndraw.draw.Cylinder(position=center,
                                                    rotation=rotation_angles,
                                                    radius_bottom=radius,
                                                    radius_top=radius,
                                                    height=length,
                                                    material=mat))

            elif shape_type == "Wall":
                dist = shape.dist
                normal = np.array(shape.normal)

                position = dist * normal
                helper = WallIntersection(
                    plane_normal=normal, plane_point=position, box_l=self.system.box_l)
                corners = helper.get_intersections()

                base_position = np.copy(corners[0])
                corners -= base_position

                # Rotate plane to align with z-axis, Custom2DShape only works
                # in the xy-plane
                unit_z = np.array([0, 0, 1])
                r, _ = scipy.spatial.transform.Rotation.align_vectors(
                    [unit_z], [normal])
                rotated_corners = r.apply(corners)

                # Sort corners in a clockwise order, except the first corner
                angles = np.arctan2(
                    rotated_corners[1:, 1], rotated_corners[1:, 0])
                sorted_indices = np.argsort(angles)
                sorted_corners = rotated_corners[1:][sorted_indices]
                sorted_corners = np.vstack(
                    [rotated_corners[0], sorted_corners])[:, :2]

                r, _ = scipy.spatial.transform.Rotation.align_vectors(
                    [normal], [unit_z])
                euler_angles = r.as_euler("xyz")

                # invert the z-axis, unsure why this is needed, maybe
                # different coordinate systems
                euler_angles[2] *= -1.

                objects.append(zndraw.draw.Custom2DShape(
                    position=base_position, rotation=euler_angles,
                    points=sorted_corners, material=mat))

            elif shape_type == "Sphere":
                center = shape.center
                radius = shape.radius

                objects.append(
                    zndraw.draw.Sphere(position=center, radius=radius, material=mat))

            elif shape_type == "Rhomboid":
                a = shape.a
                b = shape.b
                c = shape.c
                corner = shape.corner

                objects.append(
                    zndraw.draw.Rhomboid(position=corner, vectorA=a, vectorB=b, vectorC=c, material=mat))

            elif shape_type == "Ellipsoid":
                center = shape.center
                a = shape.a
                b = shape.b

                objects.append(zndraw.draw.Ellipsoid(position=center,
                               a=a, b=b, c=b, material=mat))

            else:
                raise NotImplementedError(
                    f"Shape of type {shape_type} isn't available in ZnDraw")

            self.zndraw.geometries = objects


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
