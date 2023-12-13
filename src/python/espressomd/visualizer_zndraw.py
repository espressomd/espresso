import numpy as np
from zndraw import ZnDraw
from zndraw.frame import Frame
from zndraw.vec import LatticeVecField, OriginVecField
import networkx as nx
import webbrowser
import time


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
            "vector_field_lattice_constant": None,
            "open_browser": True,
            "mozilla": False, #TODO check automatically for mozilla
            "jupyter": True
        }

        for key in kwargs:
            if key not in self.info:
                raise ValueError(f'{key} is no valid visualization property')
            else:
                self.info[key] = kwargs[key]

        self.vis = ZnDraw(jupyter=False)
        self.vis.socket.sleep(2)

        if self.info["open_browser"]:
            if self.info["mozilla"] and self.info["jupyter"]:
                # Annoying workaround for firefox not allowing new tabs to be opened by command-line
                # if there is already a window open. 
                from IPython.display import Javascript, display
                display(Javascript('window.open("{url}");'.format(url=self.vis.url)))
            else:
                webbrowser.open_new_tab(self.vis.url)
            time.sleep(4) #give the website time to open
                
    def distance(self, p1, p2):
        """
        Calculates non periodic distance between particles

        #TODO this should be optimized. 
        """
        if self.info["folded"]:
            return np.linalg.norm(p1.pos_folded - p2.pos_folded)
        else:
            return np.linalg.norm(p1.pos - p2.pos)

    def get_bonds(self):
        if not self.info["bonds"]:
            return nx.empty_graph(), False
        else:
            connectivity_matrix = np.zeros((self.length, self.length))

            for particle in self.system.part.all():
                if not particle.bonds:
                    continue
                for bond in particle.bonds:
                    bond_id = bond[-1]
                    distance = self.distance(
                        self.system.part.by_id(bond_id), particle)
                    if distance > 0.5 * min(self.system.box_l):
                        break
                    connectivity_matrix[particle.id, bond_id] = 1

            bond_graph = nx.from_numpy_array(connectivity_matrix)
            return bond_graph, True

    def get_colors(self):
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

    def get_sizes(self):
        if self.info["size"] is None:
            return np.ones((self.length))
        elif isinstance(self.info["size"], int):
            return self.info["size"] * np.ones((self.length))
        elif isinstance(self.info["size"], list) or isinstance(self.info["size"], np.ndarray):
            size_list = [""] * self.length
            try:
                for particle in self.system.part.all():
                    size_list[particle.id] = size[particle.type]
                return size_list
            except IndexError:
                raise AttributeError(
                    "The number of sizes does not match the number of particle types")

    def get_box(self):
        self.sim_box = np.array([[self.system.box_l[0], 0,0], [0,self.system.box_l[1], 0], [0,0,self.system.box_l[2]]])
        
    def get_vector_field(self):
        lbf = self.system.lb
        vectors = []
        origins = []
        lattice_constant = self.info["vector_field_lattice_constant"]
        for i in range(lattice_constant[0]):
            for j in range(lattice_constant[1]):
                for k in range(lattice_constant[2]):
                    origin = list(self.sim_box[0]/lattice_constant[0]*(i+0.5) 
                    + self.sim_box[1]/lattice_constant[1]*(j+0.5) 
                    + self.sim_box[2]/lattice_constant[2]*(k+0.5))
                    origins.append(origin)
                    vectors.append(list(lbf.get_interpolated_velocity(pos=origin)))
        vector_field = OriginVecField(vectors=vectors, origins=origins, color="#e6194B")
        return vector_field

    def system_to_frame(self):
        """
        Converts the espresso.system.System object to a zndraw.Frame object using the
        given parameters, such as particles sizes and colors.
        """

        self.length = len(self.system.part.all())

        if self.info["folded"]:
            pos = self.system.part.all().pos_folded
        else:
            pos = self.system.part.all().pos

        bond_graph, draw_bonds = self.get_bonds()
        color_list = self.get_colors()
        sizes = self.get_sizes()
        self.get_box()

        frame = Frame(positions=pos, 
                      numbers=sizes,
                      colors=color_list,
                      cell=self.sim_box,
                      connectivity=bond_graph,
                      bonds=draw_bonds,
                      auto_bonds=False,
                      vector_field=self.get_vector_field()
                      )

        return frame

    def update(self):
        new_frame = self.system_to_frame()
        self.vis[self.current_frame] = new_frame
        self.current_frame += 1
