from __future__ import print_function
import espressomd
from espressomd import thermostat
from espressomd.interactions import HarmonicBond
from espressomd.visualization_opengl import *
import numpy as np
from threading import Thread

box_l = 10
system = espressomd.System(box_l=[box_l]*3)
system.set_random_state_PRNG()
np.random.seed(seed=system.seed)

visualizer = openGLLive(system, background_color=[1, 1, 1], bond_type_radius = [0.2])

system.time_step = 0.0005
system.cell_system.skin = 0.1
#system.cell_system.min_num_cells = 1

system.box_l = [box_l, box_l, box_l]

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=2, sigma=1,
    cutoff=3, shift="auto")

system.bonded_inter[0] = HarmonicBond(k=5.0, r_0=1.0)

n_part = 200
for i in range(n_part):
    system.part.add(id=i, pos=np.random.random(3) * system.box_l)

for i in range(0,n_part - 1,2):
    system.part[i].add_bond((system.bonded_inter[0], system.part[i + 1].id))

print("E before minimization:", system.analysis.energy()["total"])
system.minimize_energy.init(f_max=0.0, gamma=30.0,
                            max_steps=10000, max_displacement=0.1)
system.minimize_energy.minimize()
print("E after minimization:", system.analysis.energy()["total"])

system.thermostat.set_npt(kT=2.0, gamma0=1.0, gammav=0.01)
system.integrator.set_isotropic_npt(ext_pressure=1.0, piston=0.01)


def main():
    p = 1.0
    b = system.box_l
    cnt = 0
    while True:
        system.integrator.run(1)
        if cnt > 1000:
            print("Pressure:", system.analysis.pressure()
                  ['total'], "Box:", system.box_l)
            cnt = 0

        visualizer.update()
        cnt += 1


# Start simulation in seperate thread
t = Thread(target=main)
t.daemon = True
t.start()

visualizer.start()
