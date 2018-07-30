#
# Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
#   Max-Planck-Institute for Polymer Research, Theory Group
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
from espressomd import assert_features, electrostatics, electrostatic_extensions
from espressomd.shapes import Wall
from espressomd.visualization_opengl import *
import numpy
from threading import Thread
from time import sleep

assert_features(["ELECTROSTATICS", "CONSTRAINTS", "MASS", "LENNARD_JONES"])

system = espressomd.System(box_l=[1.0, 1.0, 1.0])
system.seed  = system.cell_system.get_state()['n_nodes'] * [1234]
numpy.random.seed(system.seed)

print("\n--->Setup system")

# System parameters
n_part = 1000
n_ionpairs = n_part / 2
density = 1.1138
time_step = 0.001823
temp = 1198.3
gamma = 50
#l_bjerrum = 0.885^2 * e^2/(4*pi*epsilon_0*k_B*T)
l_bjerrum = 130878.0 / temp
Ez = 0

num_steps_equilibration = 500

# Particle parameters
types = {"Cl":          0, "Na": 1, "Electrode": 2}
numbers = {"Cl": n_ionpairs, "Na": n_ionpairs}
charges = {"Cl": -1.0, "Na": 1.0}
lj_sigmas = {"Cl":       3.85, "Na": 2.52,  "Electrode": 3.37}
lj_epsilons = {"Cl":     192.45, "Na": 17.44, "Electrode": 24.72}

lj_cuts = {"Cl":        3.0 * lj_sigmas["Cl"],
           "Na":        3.0 * lj_sigmas["Na"],
           "Electrode": 3.0 * lj_sigmas["Electrode"]}

masses = {"Cl":     35.453, "Na": 22.99, "Electrode": 12.01}

# Setup System
box_l = (n_ionpairs * sum(masses.values()) / density)**(1. / 3.)
box_z = box_l + 2.0 * (lj_sigmas["Electrode"])
box_volume = box_l * box_l * box_z
elc_gap = box_z * 0.15
system.box_l = [box_l, box_l, box_z + elc_gap]
system.periodicity = [1, 1, 1]
system.time_step = time_step
system.cell_system.skin = 0.3
system.thermostat.set_langevin(kT=temp, gamma=gamma)

# Visualizer
visualizer = openGLLive(system, camera_position=[-3 * box_l, box_l * 0.5, box_l * 0.5], camera_right=[
                        0, 0, 1], drag_force=5 * 298, background_color=[1, 1, 1], light_pos=[30, 30, 30], ext_force_arrows_scale=[0.0001], ext_force_arrows=False)

# Walls
system.constraints.add(shape=Wall(
    dist=0, normal=[0, 0, 1]), particle_type=types["Electrode"])
system.constraints.add(shape=Wall(
    dist=-box_z, normal=[0, 0, -1]), particle_type=types["Electrode"])

# Place particles
for i in range(int(n_ionpairs)):
    p = numpy.random.random(3) * box_l
    p[2] += lj_sigmas["Electrode"]
    system.part.add(id=len(system.part),
                    type=types["Cl"],  pos=p, q=charges["Cl"], mass=masses["Cl"])
for i in range(int(n_ionpairs)):
    p = numpy.random.random(3) * box_l
    p[2] += lj_sigmas["Electrode"]
    system.part.add(id=len(system.part),
                    type=types["Na"],  pos=p, q=charges["Na"], mass=masses["Na"])

# Lennard-Jones interactions parameters


def combination_rule_epsilon(rule, eps1, eps2):
    if rule == "Lorentz":
        return (eps1 * eps2)**0.5
    else:
        return ValueError("No combination rule defined")


def combination_rule_sigma(rule, sig1, sig2):
    if rule == "Berthelot":
        return (sig1 + sig2) * 0.5
    else:
        return ValueError("No combination rule defined")


for s in [["Cl", "Na"], ["Cl", "Cl"], ["Na", "Na"], ["Na", "Electrode"], ["Cl", "Electrode"]]:
    lj_sig = combination_rule_sigma(
        "Berthelot", lj_sigmas[s[0]], lj_sigmas[s[1]])
    lj_cut = combination_rule_sigma("Berthelot", lj_cuts[s[0]], lj_cuts[s[1]])
    lj_eps = combination_rule_epsilon(
        "Lorentz", lj_epsilons[s[0]], lj_epsilons[s[1]])

    system.non_bonded_inter[types[s[0]], types[s[1]]].lennard_jones.set_params(
        epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")

system.minimize_energy.init(
    f_max=10, gamma=10, max_steps=2000, max_displacement=0.1)
system.minimize_energy.minimize()

print("\n--->Tuning Electrostatics")
p3m = electrostatics.P3M(bjerrum_length=l_bjerrum, accuracy=1e-2)
system.actors.add(p3m)
elc = electrostatic_extensions.ELC(gap_size=elc_gap, maxPWerror=1e-3)
system.actors.add(elc)


def increaseElectricField():
    global Ez
    Ez += 1000
    for p in system.part:
        p.ext_force = [0, 0, Ez * p.q]
    print(Ez)


def decreaseElectricField():
    global Ez
    Ez -= 1000
    for p in system.part:
        p.ext_force = [0, 0, Ez * p.q]
    print(Ez)


# Register buttons
visualizer.keyboardManager.registerButton(KeyboardButtonEvent(
    'u', KeyboardFireEvent.Hold, increaseElectricField))
visualizer.keyboardManager.registerButton(KeyboardButtonEvent(
    'j', KeyboardFireEvent.Hold, decreaseElectricField))


def main():
    print("\n--->Integration")
    system.time = 0.0

    while True:
        system.integrator.run(1)
        visualizer.update()


# Start simulation in seperate thread
t = Thread(target=main)
t.daemon = True
t.start()

# Start blocking visualizer
visualizer.start()
