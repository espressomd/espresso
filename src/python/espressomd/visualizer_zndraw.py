import numpy as np
from zndraw import ZnDraw
from zndraw.frame import Frame
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

        self.vis = ZnDraw(jupyter=False)
        self.vis.socket.sleep(2)

        webbrowser.open(url=self.vis.url, new=2)
        time.sleep(1)  # give the browser time to connect before sending data

        self.info = {
            "folded": False,
            "bonds": False,
            "colors": None,
            "size": None
        }

        for key in kwargs:
            if key not in self.info:
                raise ValueError(f'{key} is no valid visualization property')
            else:
                self.info[key] = kwargs[key]

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
        return bond_graph

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

        if not self.info["bonds"]:
            draw_bonds = False
            bond_graph = nx.Graph()   
        else:
            bond_graph = self.get_bonds()
            draw_bonds = True

        if self.info["colors"] is None:
            color_list = ["#ffffff"] * self.length
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
            color_list = [self.color_dict[self.info["colors"]]] * self.length

        if self.info["size"] is None:
            size = np.ones((self.length))
        elif isinstance(self.info["size"], int):
            size = self.info["size"] * np.ones((self.length))
        elif isinstance(self.info["size"], list) or isinstance(self.info["size"], np.ndarray):
            size_list = [""] * self.length
            try:
                for particle in self.system.part.all():
                    size_list[particle.id] = size[particle.type]
                return size_list
            except IndexError:
                raise AttributeError(
                    "The number of sizes does not match the number of particle types")

        frame = Frame(positions=pos, 
                      numbers=size,
                      colors=color_list,
                      connectivity=bond_graph,
                      auto_bonds=False,
                      bonds=draw_bonds)

        return frame

    def update(self):
        present_frame = self.system_to_frame()

        self.vis[self.current_frame] = present_frame
        self.current_frame += 1
