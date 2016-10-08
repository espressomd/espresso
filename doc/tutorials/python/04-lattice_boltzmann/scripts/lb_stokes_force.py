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
from espressomd import System, lb
import numpy as np

# System setup
system = System()
agrid = 1
system.box_l = [64+2*agrid, 64+2*agrid, 64]
system.time_step = 0.1
system.cell_system.skin = 0.4

# The temperature is zero.
system.thermostat.set_lb(kT=0)

# LB Parameters
v = 0.01 # The boundary speed 
kinematic_visc = 1.0




# Invoke LB fluid
lbf = lb.LBFluid_GPU(visc=kinematic_visc, dens=1, agrid=agrid, tau=system.time_step, fric=10)
system.actors.add(lbf)

# TODO
# # Four walls make an infinite square channel along z direction
# lbboundary wall normal -1. 0. 0. dist [ expr -(1+$box_width) ] velocity 0.00 0 $v 
# lbboundary wall normal  1. 0. 0. dist 1. velocity 0. 0. $v
# lbboundary wall normal  0 -1. 0. dist [ expr -(1+$box_width) ] velocity 0.00 0 $v 
# lbboundary wall normal  0  1. 0. dist 1. velocity 0. 0. $v

radius = 5.5
# TODO
# lbboundary sphere center [ expr 0.5*($box_width+2*$agrid) ] [ expr 0.5*($box_width+2*$agrid) ] [ expr 0.5*$box_length ]\
#      radius $radius direction +1

lbf.print_vtk_boundary(boundary.vtk)

for i in range(30):
    system.integrator.run(400)
    lbf.print_vtk_velocity("fluid%04i.vtk" %i)

    # TODO
    # # print out force on the sphere
    # puts [ lindex [  lbboundary force  4 ]  2 ]

stokes_force = 6*np.pi*kinematic_visc*radius*v
print("Stokes' Law says: f=%f" %stokes_force)

