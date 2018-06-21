import espressomd
import espressomd.lb
import espressomd.shapes
import espressomd.lbboundaries
from espressomd.visualization_opengl import *

system = espressomd.System(box_l=[10.0, 10.0, 5.0])
system.time_step = 0.01
system.cell_system.skin = 0.4

lb_fluid = espressomd.lb.LBFluid(agrid=1.0, fric=1.0, dens=1.0, visc=1.0, tau=0.01, ext_force_density=[0, 0, 0.15])
system.actors.add(lb_fluid)

cylinder_shape = espressomd.shapes.Cylinder(
        center = [5.0, 5.0, 5.0],
        axis = [0, 0, 1],
        direction = -1,
        radius = 4.0,
        length = 20.0)
cylinder_boundary = espressomd.lbboundaries.LBBoundary(shape=cylinder_shape)
system.lbboundaries.add(cylinder_boundary)

visualizer = openGLLive(system, 
                        background_color = [1,1,1],
                        camera_position = [5,5,25],
                        LB_draw_boundaries = True, 
                        LB_draw_nodes = True, 
                        LB_draw_node_boundaries = True)

visualizer.run(1)
