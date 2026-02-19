#
# Copyright (C) 2025 The ESPResSo project
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
#

"""
Simulate the warm-up of a Lennard-Jones fluid by a Langevin thermostat.
An ASE calculator :cite:`hjorthlarsen17a` computes the Lennard-Jones force
as a time-dependent external force
(particle property :attr:`~espressomd.particle_data.ParticleHandle.ext_force`)
through the :class:`ASE interface <espressomd.plugins.ase.ASEInterface>`.
"""

import numpy as np
import espressomd
import espressomd.plugins.ase
from ase.calculators.lj import LennardJones

required_features = ["LENNARD_JONES", "EXTERNAL_FORCES"]
espressomd.assert_features(required_features)

# System parameters
box_l = 10.7437
density = 0.7
system = espressomd.System(box_l=[box_l] * 3)
np.random.seed(seed=42)

# Interaction parameters (purely repulsive Lennard-Jones)
lj_eps = 1.
lj_sig = 1.
lj_cut = 2.**(1 / 6) * lj_sig

# Integration parameters
system.time_step = 0.005
system.cell_system.skin = 0.4

minimize_steps = 10
int_steps = 40
int_n_times = 16

# Particle setup
volume = system.volume()
n_part = int(volume * density)

system.part.add(pos=np.random.random((n_part, 3)) * system.box_l)

## Setup of ASE interface and Lennard-Jones calculator

print("ASE Lennard-Jones calculator setup")

ase = espressomd.plugins.ase.ASEInterface(system, system.part.all())

# ASE calculator tor provide Lennard-Jones forces
lj = LennardJones(sigma=lj_sig, epsilon=lj_eps, rc=lj_cut, smooth=False)
ase.atoms.calc = lj

# Overlap removal via steepest descent
print("Overlap removal")
system.integrator.set_steepest_descent(f_max=0., gamma=1e-3,
                                       max_displacement=lj_sig / 100.)
while True:
    ase.integrate(minimize_steps, lj)
    max_force = np.amax(np.linalg.norm(system.part.all().f, axis=1))
    print(f"{max_force=:.2e}")
    if max_force < 10.:
        break

# Normal integration
system.integrator.set_symplectic_euler()
system.integrator.run(0, reuse_forces=False)
forces = system.part.all().f
ase.update_ase()
box_l = np.copy(system.box_l)
ase_pos = np.remainder(np.array([a.position for a in ase.atoms]), box_l)
np.testing.assert_allclose(
    np.remainder(ase_pos - np.copy(system.part.all().pos), box_l), 0.)

print("LJ fluid warm-up to T=1 (in MD units)")
system.thermostat.set_langevin(kT=1., gamma=1., seed=42)
instantaneous_temperatures = []
for i in range(int_n_times):
    ase.integrate(int_steps, lj)
    energy = system.analysis.energy()
    T_inst = 2. / 3. * energy["kinetic"] / n_part
    instantaneous_temperatures.append(T_inst)
    print(f"{T_inst=:.2f}")
