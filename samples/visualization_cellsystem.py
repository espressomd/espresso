# Copyright (C) 2010-2019 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
Visualize the system cells and MPI domains. Run ESPResSo in parallel
to color particles by node. With OpenMPI, this can be achieved using
``mpiexec -n 4 ./pypresso ../samples/visualization_cellsystem.py``.
Set property ``system.cell_system.node_grid = [i, j, k]`` (with ``i * j * k``
equal to the number of MPI ranks) to change the way the cellsystem is
partitioned. Only the domain of MPI rank 0 will be shown in wireframe.
"""

import espressomd
from espressomd.minimize_energy import steepest_descent
from espressomd.visualization_opengl import openGLLive
import numpy as np

required_features = ["LENNARD_JONES"]
espressomd.assert_features(required_features)

box = [40, 30, 20]
system = espressomd.System(box_l=box)
visualizer = openGLLive(
    system,
    window_size=[800, 800],
    background_color=[0, 0, 0],
    camera_position=[20, 15, 80],
    particle_coloring='node',
    draw_nodes=True,
    draw_cells=True)

system.time_step = 0.0005
system.cell_system.set_domain_decomposition(use_verlet_lists=True)
system.cell_system.skin = 0.4
#system.cell_system.node_grid = [i, j, k]

for i in range(100):
    system.part.add(pos=box * np.random.random(3))

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=100.0, sigma=1.0, cutoff=3.0, shift="auto")

energy = system.analysis.energy()
print("Before Minimization: E_total = {}".format(energy['total']))
steepest_descent(system, f_max=50, gamma=30.0, max_steps=10000,
                 max_displacement=0.001)
energy = system.analysis.energy()
print("After Minimization: E_total = {}".format(energy['total']))

print("Tune skin")
system.cell_system.tune_skin(0.1, 4.0, 1e-1, 1000)
print(system.cell_system.get_state())

system.thermostat.set_langevin(kT=1, gamma=1, seed=42)

visualizer.run(1)
