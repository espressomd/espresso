from __future__ import print_function
import espressomd
from espressomd import thermostat
from espressomd import integrate
from espressomd.interactions import HarmonicBond
from espressomd import visualization
import numpy
from threading import Thread

box_l = 50
n_part = 200

system = espressomd.System()
system.time_step = 0.01
system.cell_system.skin = 0.4
system.thermostat.set_langevin(kT=0.1, gamma=20.0)

system.box_l = [box_l, box_l, box_l]

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=0, sigma=1,
    cutoff=2, shift="auto")
system.bonded_inter[0] = HarmonicBond(k=0.5, r_0=1.0)
system.bonded_inter[1] = HarmonicBond(k=0.5, r_0=1.0)

for i in range(n_part):
    system.part.add(id=i, pos=numpy.random.random(3) * system.box_l)

for i in range(n_part - 1):
    system.part[i].add_bond((system.bonded_inter[0], system.part[i + 1].id))

#visualizer = visualization.mayaviLive(system)
visualizer = visualization.openGLLive(system, bond_type_radius=[0.3])

system.minimize_energy.init(
    f_max=10, gamma=50.0, max_steps=1000, max_displacement=0.2)
system.minimize_energy.minimize()

def main():
    while True:
        system.integrator.run(1)
        visualizer.update()


# Start simulation in seperate thread
t = Thread(target=main)
t.daemon = True
t.start()

visualizer.start()
