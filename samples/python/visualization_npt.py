from __future__ import print_function
import espressomd
from espressomd import thermostat
from espressomd.visualization_opengl import *
import numpy
from threading import Thread

system=espressomd.System()
visualizer = openGLLive(system, background_color=[1,1,1])

system.time_step = 0.005
system.cell_system.skin = 0.4

box_l=10
system.box_l = [box_l, box_l, box_l]

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=1, sigma=1,
    cutoff=2**(1./6.), shift="auto")

for i in range(100):
    system.part.add(id=i, pos=numpy.random.random(3) * system.box_l)

print("E before minimization:",system.analysis.energy()["total"])
system.minimize_energy.init(f_max = 0.0, gamma = 30.0, max_steps = 10000, max_displacement= 0.1)
system.minimize_energy.minimize()
print("E after minimization:",system.analysis.energy(i)["total"])

system.thermostat.set_npt(kT=1.0, gamma0=1.0, gammav=0.1)
system.integrator.set_isotropic_npt(ext_pressure=1.0, piston=1.0)

def main():
    p=1.0
    b=system.box_l
    while True:
        system.integrator.run(1)
        print("Pressure:", system.analysis.pressure()['total'], "Box:", system.box_l)
        visualizer.update()

#Start simulation in seperate thread
t = Thread(target=main)
t.daemon = True
t.start()

visualizer.start()
