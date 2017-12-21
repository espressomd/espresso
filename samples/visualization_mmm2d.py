from __future__ import print_function
from espressomd import *
from espressomd.shapes import *
from espressomd import electrostatics
import numpy
from threading import Thread
from math import *
from espressomd.visualization_opengl import *

system = espressomd.System()
visualizer = openGLLive(system, constraint_type_colors= [[1,1,1,1]], camera_position = [50,15,15], camera_right = [0,0,-1] )


box_l = 20
system.box_l = [box_l, box_l, box_l]
system.time_step = 0.02
system.cell_system.skin = 0.4
system.cell_system.set_layered(n_layers=5,use_verlet_lists = False)
system.periodicity = [1, 1, 0]

qion = 1
for i in range(300):
    rpos = numpy.random.random(3) * box_l
    system.part.add(pos=rpos, type = 0, q = qion)
    qion *= -1

system.constraints.add(shape=Wall(dist=0,normal=[0,0,1]),particle_type=1)
system.constraints.add(shape=Wall(dist=-box_l,normal=[0,0,-1]),particle_type=1)
        
WCA_cut = 2.**(1. / 6.)
system.non_bonded_inter[0,1].lennard_jones.set_params(epsilon=1.0, sigma=1.0, cutoff=WCA_cut, shift="auto")
system.non_bonded_inter[0,0].lennard_jones.set_params(epsilon=1.0, sigma=1.0, cutoff=WCA_cut, shift="auto")

energy = system.analysis.energy()
print("Before Minimization: E_total=", energy['total'])
system.minimize_energy.init(f_max = 10, gamma = 50.0, max_steps = 1000, max_displacement= 0.2)
system.minimize_energy.minimize()
energy = system.analysis.energy()
print("After Minimization: E_total=", energy['total'])

system.thermostat.set_langevin(kT=0.1, gamma=1.0)

mmm2d = electrostatics.MMM2D(prefactor = 10.0, maxPWerror = 1e-3, const_pot = 1, pot_diff = 50.0)
system.actors.add(mmm2d)

def main():

    while True:
        system.integrator.run(1)
        visualizer.update()


# Start simulation in seperate thread
t = Thread(target=main)
t.daemon = True
t.start()

# Start blocking visualizer
visualizer.start()
