import numpy as np
from zndraw import ZnDraw
from znframe.frame import Frame
from znframe.vec import LatticeVecField, OriginVecField
from zndraw.utils import SphereGeometry, CylinderGeometry, PlaneGeometry
import networkx as nx
import webbrowser
import time
import espressomd.shapes


class Visualizer():

    color_dict = {"black": "#000000",
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
                  "white": "#ffffff"}

    def __init__(self, system, **kwargs):
        self.system = system
        self.current_frame = 0

        self.info = {
            "folded": False,
            "bonds": False,
            "colors": None,
            "size": None,
            "vector_field": None,
            "open_browser": True,
            "jupyter": True
        }

        for key in kwargs:
            if key not in self.info:
                raise ValueError(f'{key} is no valid visualization property')
            else:
                self.info[key] = kwargs[key]

        self.vis = ZnDraw(jupyter=False)
        self.vis.socket.sleep(2)

        if self.info["jupyter"]:
            from IPython.display import display, IFrame
            display(IFrame(src=self.vis.url, width="100%", height="700px"))
            time.sleep(4)
        elif self.info["open_browser"]:
            if webbrowser.get().name == "firefox":
                print("browser is:", webbrowser.get().name)
                from IPython.display import Javascript, display
                display(
                    Javascript(
                        'window.open("{url}");'.format(
                            url=self.vis.url)))
            else:
                webbrowser.open_new_tab(self.vis.url)
            time.sleep(4)  # give the website time to open

        self._get_box()

    def _distance(self, p1, p2):
        """
        Calculates non periodic distance between particles
        """
        if self.info["folded"]:
            return np.linalg.norm(p1.pos_folded - p2.pos_folded)
        else:
            return np.linalg.norm(p1.pos - p2.pos)

    def _get_bonds(self):
        if not self.info["bonds"]:
            return []
        else:
            bonds = []

            for particle in self.system.part.all():
                if not particle.bonds:
                    continue
                for bond in particle.bonds:
                    bond_id = bond[-1]
                    distance = self._distance(
                        self.system.part.by_id(bond_id), particle)
                    if distance > 0.5 * min(self.system.box_l):
                        break
                    bonds.append((particle.id, bond_id, 1))
            return bonds

    def _get_colors(self):
        if self.info["colors"] is None:
            return ["#ffffff"] * self.length
        elif isinstance(self.info["colors"], list):
            color_list = [""] * self.length
            try:
                for particle in self.system.part.all():
                    color_list[particle.id] = self.color_dict[self.info["colors"]
                                                              [particle.type]]
            except IndexError:
                raise AttributeError(
                    "The number of colors does not match the number of particle types")
        elif isinstance(self.info["colors"], str):
            return [self.color_dict[self.info["colors"]]] * self.length

    def _get_sizes(self):
        if self.info["size"] is None:
            return np.ones((self.length))
        elif isinstance(self.info["size"], int):
            return self.info["size"] * np.ones((self.length), dtype=int)
        elif isinstance(self.info["size"], list) or isinstance(self.info["size"], np.ndarray):
            size_list = [""] * self.length
            try:
                for particle in self.system.part.all():
                    size_list[particle.id] = int(size[particle.type])
                return size_list
            except IndexError:
                raise AttributeError(
                    "The number of sizes does not match the number of particle types")

    def _get_box(self):
        self.sim_box = np.array([[self.system.box_l[0], 0, 0], [
                                0, self.system.box_l[1], 0], [0, 0, self.system.box_l[2]]])

    def draw_constraints(self):
        self.shapes = []

        for constraint in self.system.constraints:
            if isinstance(constraint,
                          espressomd.constraints.ShapeBasedConstraint):
                shape = constraint.get_parameter("shape")
                shape_name = shape.name()

                if shape_name == "Shapes::Cylinder":
                    center = shape.center
                    axis = shape.axis
                    length = shape.length
                    radius = shape.radius

                    self.vis.draw_obj(
                        CylinderGeometry(
                            position=list(center),
                            direction=list(axis),
                            height=length,
                            radius=radius,
                            color="#000000",
                            opacity=0.2,
                            locked=True,
                            wireframe=True))

                elif shape_name == "Shapes::Wall":
                    dist = shape.dist
                    normal = np.array(shape.normal)

                    self.vis.draw_obj(PlaneGeometry(position=list(dist * normal), direction=list(normal), color="#000000",
                                                    width=np.max(np.diag(self.sim_box)), height=np.max(np.diag(self.sim_box)), opacity=0.2, locked=True, wireframe=True))

                elif shape_name == "Shapes::Sphere":
                    center = shape.center
                    radius = shape.radius

                    self.vis.draw_obj(
                        SphereGeometry(
                            position=list(center),
                            radius=radius,
                            color="#000000",
                            opacity=0.2,
                            locked=True,
                            wireframe=True))

    def _system_to_frame(self):
        """
        Converts the espresso.system.System object to a zndraw.Frame object using the
        given parameters, such as particles sizes and colors.
        """

        self.length = len(self.system.part.all())

        if self.info["folded"]:
            pos = self.system.part.all().pos_folded
        else:
            pos = self.system.part.all().pos

        bond_graph = self._get_bonds()
        color_list = self._get_colors()
        sizes = self._get_sizes()

        if self.info["vector_field"] is not None:
            vec = self.info["vector_field"]()
        else:
            vec = []

        arrays = {"colors": color_list}

        frame = Frame(positions=pos, 
                      numbers=sizes,
                      arrays=arrays,
                      cell=self.sim_box,
                      connectivity=bond_graph,
                      vector_field=vec,
                      recompute=[]
                      )
        return frame

    def update(self):

        new_frame = self._system_to_frame()
        self.vis[self.current_frame] = new_frame
        self.current_frame += 1


class LBField:

    def __init__(self, system, lattice_constants, scale=1):
        self.system = system
        self.lc = lattice_constants
        self.scale = scale

    def __call__(self):
        lbf = self.system.lb
        sim_box = np.array([[self.system.box_l[0], 0, 0], [
                           0, self.system.box_l[1], 0], [0, 0, self.system.box_l[2]]])
        vectors = []
        origins = []

        for i in range(self.lc[0]):
            for j in range(self.lc[1]):
                for k in range(self.lc[2]):
                    origin = list(sim_box[0] / (self.lc[0] + 1) * (i + 1) 
                                  + sim_box[1] / (self.lc[1] + 1) * (j + 1) 
                                  + sim_box[2] / (self.lc[2] + 1) * (k + 1))
                    origins.append(origin)
                    if self.scale != 1:
                        vectors.append(
                            list(
                                self.scale *
                                np.array(
                                    lbf.get_interpolated_velocity(
                                        pos=origin))))
                    else:
                        vectors.append(
                            list(
                                lbf.get_interpolated_velocity(
                                    pos=origin)))

        return OriginVecField(
            vectors=vectors, origins=origins, color="#e6194B")


# typing should be addressed, as switching between numpy and list is not
# optimal
class LatticeVectorField:

    def __init__(self, vectors, origin, color,
                 box, lattice_constants, scale=1):
        self.vectors = vectors
        self.origin = origin
        self.color = color
        self.box = box
        self.lc = lattice_constants
        self.scale = scale

    def __call__(self):
        vectors = self.vectors
        if self.scale != 1:
            vectors = self.scale * np.array(vectors)
        return LatticeVecField(vectors=list(
            self.vectors), origin=self.origin, color=self.color, box=self.box, density=self.lc)


class OriginVectorField:

    def __init__(self, vectors, origins, scale=1):
        self.vectors = vectors
        self.origins = origins
        self.scale = scale

    def __call__(self):
        vectors = self.vectors
        if self.scale != 1:
            vectors = self.scale * np.array(vectors)
        return OriginVecField(vectors=list(self.vectors), origins=self.origins)
