### Measuring the force on a single sphere immersed in a fluid with
### fixed velocity boundary conditions created by two 
### walls at finite distance.
### The force is compared to th analytical result F=6 pi eta r v
### i.e. the stokes force on the particles.


# We create a box of size box_width x box_width x box_length and 
# place an object in the center. We measure the drag force 
# in z direction. We create walls in the xz and yz plane at the box 
# boundaries, where the velocity is fixed to $v. 
#
from espressomd import lb, lbboundaries, shapes
import numpy as np
import sys
# System setup
system = espressomd.System(box_l=[1.0, 1.0, 1.0])
system.seed  = system.cell_system.get_state()['n_nodes'] * [1234]
agrid = 1
radius = 5.5
box_width = 64
real_width = box_width+2*agrid
box_length = 64
system.box_l = [real_width, real_width, box_length]
system.time_step = 0.2
system.cell_system.skin = 0.4

# The temperature is zero.
system.thermostat.set_lb(kT=0)

# LB Parameters
v = [0,0,0.01] # The boundary slip 
kinematic_visc = 1.0

# Invoke LB fluid
lbf = lb.LBFluid(visc=kinematic_visc, dens=1, agrid=agrid, tau=system.time_step, fric=1)
system.actors.add(lbf)

# Setup walls
walls = [None] * 4
walls[0] = lbboundaries.LBBoundary(shape=shapes.Wall(normal=[-1,0,0], 
                                   dist = -(1+box_width)), velocity=v)
walls[1] = lbboundaries.LBBoundary(shape=shapes.Wall(normal=[1,0,0], dist = 1),
                                 velocity=v)
walls[2] = lbboundaries.LBBoundary(shape=shapes.Wall(normal=[0,-1,0], 
                                   dist = -(1+box_width)), velocity=v)
walls[3] = lbboundaries.LBBoundary(shape=shapes.Wall(normal=[0,1,0], dist = 1),
                                   velocity=v)

for wall in walls:
    system.lbboundaries.add(wall)

# setup sphere without slip in the middle
sphere = lbboundaries.LBBoundary(shape=shapes.Sphere(radius=radius,
                                 center = [real_width/2] * 2 + [box_length/2],
                                 direction = 1))

system.lbboundaries.add(sphere)


lbf.print_vtk_boundary("./boundary.vtk")


def size(vector):
    tmp = 0
    for k in vector:
        tmp+=k*k
    return np.sqrt(tmp)

for i in range(40):
    system.integrator.run(400)
    lbf.print_vtk_velocity("fluid%04i.vtk" %i)
    sys.stdout.write("\rloop: %02i"%i)
    sys.stdout.flush()


print("\nIntegration finished.")

# get force that is exerted on the sphere
force = sphere.get_force()
print("Measured force: f=%f" %size(force))

stokes_force = 6*np.pi*kinematic_visc*radius*size(v)
print("Stokes' Law says: f=%f" %stokes_force)

