## A script to simulate planar Poisseuille flow in Espresso
from __future__ import print_function
from espressomd import System, lb, shapes, lbboundaries
import numpy as np
from threading import Thread
from espressomd.visualization_opengl import *

# System setup
system = System()
system.time_step = 0.01
system.cell_system.skin = 0.2

box_l = 16
system.box_l = [box_l] * 3

visualizer = openGLLive(system, 
                        LB = True, 
                        LB_plane_dist = 8, 
                        LB_plane_axis = 1, 
                        LB_vel_scale = 1e2, 
                        LB_plane_ngrid = 15, 
                        camera_position = [8,16,50], 
                        velocity_arrows=True, 
                        velocity_arrows_scale=[20.], 
                        arrows_radius= 0.1, 
                        velocity_arrows_type_colors=[[0,1,0,0.5]] )

lbf = lb.LBFluid(agrid=1.0, fric = 1.0, dens=1.0, visc=1.0, tau=0.1, ext_force=[0, 0.003, 0])
system.actors.add(lbf)
system.thermostat.set_lb(kT=0)

# Setup boundaries
walls = [lbboundaries.LBBoundary() for k in range(2)]
walls[0].set_params(shape=shapes.Wall(normal=[1,0,0], dist = 1.5))
walls[1].set_params(shape=shapes.Wall(normal=[-1,0,0], dist = -14.5))

for i in range(100):
    system.part.add(pos=np.random.random(3)*system.box_l)

for wall in walls:
    system.lbboundaries.add(wall)

def main():
    while True:
        system.integrator.run(1)
        visualizer.update()

# Start simulation in seperate thread
t = Thread(target=main)
t.daemon = True
t.start()

visualizer.start()
