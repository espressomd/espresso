""" Visualization sample for Poisseuille flow with Lattice Boltzmann.
"""

from __future__ import print_function
from espressomd import System, lb, shapes, lbboundaries
import numpy as np
from threading import Thread
from espressomd.visualization_opengl import *

required_features = ["LB", "LB_BOUNDARIES", "EXTERNAL_FORCES"]
espressomd.assert_features(required_features)

# System setup
box_l = 16
system = System(box_l=[box_l] * 3)
system.set_random_state_PRNG()
np.random.seed(seed=system.seed)

system.time_step = 0.01
system.cell_system.skin = 0.2

visualizer = openGLLive(system,
                        LB_draw_boundaries=True,
                        LB_draw_velocity_plane=True,
                        LB_plane_dist=8,
                        LB_plane_axis=1,
                        LB_vel_scale=1e2,
                        LB_plane_ngrid=15,
                        camera_position=[8, 16, 50],
                        velocity_arrows=True,
                        velocity_arrows_type_scale=[20.],
                        velocity_arrows_type_radii=[0.1],
                        velocity_arrows_type_colors=[[0, 1, 0]])

lbf = lb.LBFluid(agrid=1.0, fric=1.0, dens=1.0,
                 visc=1.0, tau=0.1, ext_force_density=[0, 0.003, 0])
system.actors.add(lbf)
system.thermostat.set_lb(kT=0)

# Setup boundaries
walls = [lbboundaries.LBBoundary() for k in range(2)]
walls[0].set_params(shape=shapes.Wall(normal=[1, 0, 0], dist=1.5))
walls[1].set_params(shape=shapes.Wall(normal=[-1, 0, 0], dist=-14.5))

for i in range(100):
    system.part.add(pos=np.random.random(3) * system.box_l)

for wall in walls:
    system.lbboundaries.add(wall)

visualizer.run(1)
