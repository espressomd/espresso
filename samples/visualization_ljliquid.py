#
# Copyright (C) 2013-2019 The ESPResSo project
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
Visualize a Lennard-Jones liquid with live plotting via matplotlib.
"""

import numpy as np
import matplotlib.pyplot as plt
from threading import Thread
import espressomd
import espressomd.visualization

required_features = ["LENNARD_JONES"]
espressomd.assert_features(required_features)

print("""
=======================================================
=                    lj_liquid.py                     =
=======================================================
""")

# System parameters
#############################################################

# 10 000  Particles
box_l = 10.7437
density = 0.7

# Interaction parameters (repulsive Lennard-Jones)
#############################################################

lj_eps = 1.0
lj_sig = 1.0
lj_cut = 1.12246

# Integration parameters
#############################################################
system = espressomd.System(box_l=[box_l] * 3)
np.random.seed(seed=42)

system.time_step = 0.001
system.cell_system.skin = 0.4

# warmup integration (steepest descent)
warm_steps = 20
warm_n_times = 30
# convergence criterion (particles are separated by at least 90% sigma)
min_dist = 0.9 * lj_sig

# integration
int_steps = 10
int_n_times = 50000


#############################################################
#  Setup System                                             #
#############################################################

# Interaction setup
#############################################################
system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")

print("LJ-parameters:")
print(system.non_bonded_inter[0, 0].lennard_jones.get_params())

# Particle setup
#############################################################

volume = box_l**3
n_part = int(volume * density)

for i in range(n_part):
    system.part.add(pos=np.random.random(3) * system.box_l)

print(
    f"Simulate {n_part} particles in a cubic box {box_l} at density {density}.")
print("Interactions:\n")
act_min_dist = system.analysis.min_dist()
print(f"Start with minimal distance {act_min_dist}")

visualizer = espressomd.visualization.openGLLive(system)

#############################################################
#  Warmup Integration                                       #
#############################################################

print(f"""\
Start warmup integration:
At maximum {warm_n_times} times {warm_steps} steps
Stop if minimal distance is larger than {min_dist}""")
print(system.non_bonded_inter[0, 0].lennard_jones)

# minimize energy using min_dist as the convergence criterion
system.integrator.set_steepest_descent(f_max=0, gamma=1e-3,
                                       max_displacement=lj_sig / 100)
i = 0
while i < warm_n_times and system.analysis.min_dist() < min_dist:
    print(f"minimization: {system.analysis.energy()['total']:+.2e}")
    system.integrator.run(warm_steps)
    i += 1
    visualizer.update()

print(f"minimization: {system.analysis.energy()['total']:+.2e}")
print()
system.integrator.set_vv()

# activate thermostat
system.thermostat.set_langevin(kT=1.0, gamma=1.0, seed=42)

#############################################################
#      Integration                                          #
#############################################################
print("\nStart integration: run %d times %d steps" % (int_n_times, int_steps))

# print initial energies
energies = system.analysis.energy()
print(energies)

plot, = plt.plot([0], [energies['total']], label="total")
plt.xlabel("Time")
plt.ylabel("Energy")
plt.legend()
plt.show(block=False)


def main_loop():
    global energies
    print(f"run at time={system.time:.2f}")

    system.integrator.run(int_steps)
    visualizer.update()

    energies = system.analysis.energy()
    plot.set_xdata(np.append(plot.get_xdata(), system.time))
    plot.set_ydata(np.append(plot.get_ydata(), energies['total']))


def main_thread():
    for _ in range(int_n_times):
        main_loop()


last_plotted = 0


def update_plot():
    global last_plotted
    current_time = plot.get_xdata()[-1]
    if last_plotted == current_time:
        return
    last_plotted = current_time
    plt.xlim(0, plot.get_xdata()[-1])
    plt.ylim(plot.get_ydata().min(), plot.get_ydata().max())
    plt.draw()
    plt.pause(0.01)


t = Thread(target=main_thread)
t.daemon = True
t.start()
visualizer.register_callback(update_plot, interval=1000)
visualizer.start()

# terminate program
print("\nFinished.")
