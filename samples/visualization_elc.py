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
Visualize charged particles confined between two plates of a capacitor with
a potential difference. The system is periodic in the *xy*-plane but has a gap
in the *z*-direction. The ELC method subtracts the electrostatic contribution
from the periodic images in the *z*-direction. The system total charge is zero.
For more details, see :ref:`Electrostatic Layer Correction (ELC)`.
"""

import numpy as np

import espressomd
import espressomd.shapes
from espressomd.minimize_energy import steepest_descent
from espressomd import electrostatics
from espressomd import electrostatic_extensions
from espressomd import visualization

required_features = ["P3M", "WCA"]
espressomd.assert_features(required_features)

box_l = 20
elc_gap = 10
potential_diff = -3.
system = espressomd.System(box_l=[box_l, box_l, box_l + elc_gap])
np.random.seed(seed=42)
visualizer = visualization.openGLLive(
    system,
    background_color=[1, 1, 1],
    constraint_type_colors=[[0, 0, 0]],
    camera_position=[70, 10, 15],
    camera_right=[0, 0, -1])

system.time_step = 0.02
system.cell_system.skin = 0.4
system.periodicity = [1, 1, 1]

qion = 1
for i in range(300):
    rpos = np.random.random(3) * box_l
    system.part.add(pos=rpos, type=0, q=qion)
    qion *= -1

system.constraints.add(shape=espressomd.shapes.Wall(
    dist=0, normal=[0, 0, 1]), particle_type=1)
system.constraints.add(shape=espressomd.shapes.Wall(
    dist=-box_l, normal=[0, 0, -1]), particle_type=1)

system.non_bonded_inter[0, 1].wca.set_params(epsilon=1.0, sigma=1.0)
system.non_bonded_inter[0, 0].wca.set_params(epsilon=1.0, sigma=1.0)

energy = system.analysis.energy()
print("Before Minimization: E_total=", energy['total'])
steepest_descent(system, f_max=10, gamma=50.0, max_steps=1000,
                 max_displacement=0.2)
energy = system.analysis.energy()
print("After Minimization: E_total=", energy['total'])

system.thermostat.set_langevin(kT=0.1, gamma=1.0, seed=42)

p3m = electrostatics.P3M(prefactor=1.0, accuracy=1e-2)

system.actors.add(p3m)

elc = electrostatic_extensions.ELC(maxPWerror=1.0, gap_size=elc_gap,
                                   const_pot=True, pot_diff=potential_diff)
system.actors.add(elc)

visualizer.run(1)
