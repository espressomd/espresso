""" Visualization of the system cells and MPI domains.
"""

import espressomd
from espressomd.visualization_opengl import openGLLive
import numpy as np

required_features = ["LENNARD_JONES"]
espressomd.assert_features(required_features)

box = [40, 30, 20]
system = espressomd.System(box_l=box)
system.set_random_state_PRNG()
visualizer = openGLLive(system, window_size=[800, 800], background_color=[0, 0, 0], camera_position=[
                        20, 15, 80], particle_coloring='node', draw_nodes=True, draw_cells=True)

system.time_step = 0.0005
system.cell_system.set_domain_decomposition(use_verlet_lists=True)
system.cell_system.skin = 0.4

for i in range(100):
    system.part.add(pos=box * np.random.random(3))

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=100.0, sigma=1.0,
    cutoff=3.0, shift="auto")

energy = system.analysis.energy()
print("Before Minimization: E_total = {}".format(energy['total']))
system.minimize_energy.init(
    f_max=50, gamma=30.0, max_steps=10000, max_displacement=0.001)
system.minimize_energy.minimize()
energy = system.analysis.energy()
print("After Minimization: E_total = {}".format(energy['total']))

print("Tune skin")
system.cell_system.tune_skin(0.1, 4.0, 1e-1, 1000)
print(system.cell_system.get_state())

system.thermostat.set_langevin(kT=1, gamma=1)

visualizer.run(1)
