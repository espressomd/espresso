#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
from __future__ import print_function
import numpy
from matplotlib import pyplot
from threading import Thread
import espressomd
from espressomd import thermostat
from espressomd import integrate
from espressomd import visualization

print("""
=======================================================
=                    lj_liquid.py                     =
=======================================================

Program Information:""")
print(espressomd.features())

dev = "cpu"

# System parameters
#############################################################

# 10 000  Particles
box_l = 10.7437
density = 0.7

# Interaction parameters (repulsive Lennard Jones)
#############################################################

lj_eps = 1.0
lj_sig = 1.0
lj_cut = 1.12246
lj_cap = 20

# Integration parameters
#############################################################
system = espressomd.System()
system.time_step = 0.001
system.cell_system.skin = 0.4
#es._espressoHandle.Tcl_Eval('thermostat langevin 1.0 1.0')
system.thermostat.set_langevin(kT=1.0, gamma=1.0)

# warmup integration (with capped LJ potential)
warm_steps = 100
warm_n_times = 30
# do the warmup until the particles have at least the distance min__dist
min_dist = 0.9

# integration
int_steps = 10
int_n_times = 50000


#############################################################
#  Setup System                                             #
#############################################################

# Interaction setup
#############################################################

system.box_l = [box_l, box_l, box_l]

system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig,
    cutoff=lj_cut, shift="auto")
system.force_cap = lj_cap

print("LJ-parameters:")
print(system.non_bonded_inter[0, 0].lennard_jones.get_params())

# Particle setup
#############################################################

volume = box_l * box_l * box_l
n_part = int(volume * density)

for i in range(n_part):
    system.part.add(id=i, pos=numpy.random.random(3) * system.box_l)

system.analysis.distto(0)

print("Simulate {} particles in a cubic simulation box {} at density {}."
      .format(n_part, box_l, density).strip())
print("Interactions:\n")
act_min_dist = system.analysis.mindist()
print("Start with minimal distance {}".format(act_min_dist))

system.cell_system.max_num_cells = 2744

# Switch between openGl/Mayavi
#visualizer = visualization.mayaviLive(system)
visualizer = visualization.openGLLive(system)

#############################################################
#  Warmup Integration                                       #
#############################################################

print("""
Start warmup integration:
At maximum {} times {} steps
Stop if minimal distance is larger than {}
""".strip().format(warm_n_times, warm_steps, min_dist))

# set LJ cap
lj_cap = 20
system.force_cap = lj_cap
print(system.non_bonded_inter[0, 0].lennard_jones)

# Warmup Integration Loop
i = 0
while (i < warm_n_times and act_min_dist < min_dist):
    system.integrator.run(warm_steps)
    # Warmup criterion
    act_min_dist = system.analysis.mindist()
#  print("\rrun %d at time=%f (LJ cap=%f) min dist = %f\r" % (i,system.time,lj_cap,act_min_dist), end=' ')
    i += 1

#   Increase LJ cap
    lj_cap = lj_cap + 10
    system.force_cap = lj_cap
    visualizer.update()

# Just to see what else we may get from the c code
# print("""
# ro variables:
# cell_grid     {0.cell_grid}
# cell_size     {0.cell_size}
# local_box_l   {0.local_box_l}
# max_cut       {0.max_cut}
# max_part      {0.max_part}
# max_range     {0.max_range}
# max_skin      {0.max_skin}
# n_nodes       {0.n_nodes}
# n_part        {0.n_part}
# n_part_types  {0.n_part_types}
# periodicity   {0.periodicity}
# verlet_reuse  {0.verlet_reuse}
#""".format(system))

#############################################################
#      Integration                                          #
#############################################################
print("\nStart integration: run %d times %d steps" % (int_n_times, int_steps))

# remove force capping
lj_cap = 0
system.force_cap = lj_cap
print(system.non_bonded_inter[0, 0].lennard_jones)

# print initial energies
energies = system.analysis.energy()
print(energies)

plot, = pyplot.plot([0], [energies['total']], label="total")
pyplot.xlabel("Time")
pyplot.ylabel("Energy")
pyplot.legend()
pyplot.show(block=False)

j = 0


def main_loop():
    global energies
    print("run %d at time=%f " % (i, system.time))

    system.integrator.run(int_steps)
    visualizer.update()

    energies = system.analysis.energy()
    plot.set_xdata(numpy.append(plot.get_xdata(), system.time))
    plot.set_ydata(numpy.append(plot.get_ydata(), energies['total']))

def main_thread():
    for i in range(0, int_n_times):
        main_loop()


last_plotted = 0


def update_plot():
    global last_plotted
    current_time = plot.get_xdata()[-1]
    if last_plotted == current_time:
        return
    last_plotted = current_time
    pyplot.xlim(0, plot.get_xdata()[-1])
    pyplot.ylim(plot.get_ydata().min(), plot.get_ydata().max())
    pyplot.draw()
    pyplot.pause(0.01)


t = Thread(target=main_thread)
t.daemon = True
t.start()
visualizer.registerCallback(update_plot, interval=1000)
visualizer.start()

# terminate program
print("\nFinished.")
