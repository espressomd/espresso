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
from espressomd.visualization_opengl import *
import numpy
from threading import Thread
from time import sleep

assert_features(["ELECTROSTATICS", "MASS", "LENNARD_JONES"])

system = espressomd.System(box_l=[1.0, 1.0, 1.0])
system.seed  = system.cell_system.get_state()['n_nodes'] * [1234]
numpy.random.seed(system.seed)

visualizer = openGLLive(system, drag_force=5 * 298,
                        background_color=[1, 1, 1], light_pos=[30, 30, 30])

# Callbacks to control temperature
temperature = 298.0


def increaseTemp():
    global temperature
    temperature += 10
    system.thermostat.set_langevin(kT=temperature, gamma=1.0)
    print(temperature)


def decreaseTemp():
    global temperature
    temperature -= 10

    if temperature > 0:
        system.thermostat.set_langevin(kT=temperature, gamma=1.0)
    else:
        temperature = 0
        system.thermostat.turn_off()
    print(temperature)


# Register buttons
visualizer.keyboardManager.registerButton(
    KeyboardButtonEvent('t', KeyboardFireEvent.Hold, increaseTemp))
visualizer.keyboardManager.registerButton(
    KeyboardButtonEvent('g', KeyboardFireEvent.Hold, decreaseTemp))


def main():

    print("\n--->Setup system")

# System parameters
    n_ppside = 10
    n_part = int(n_ppside**3)
    n_ionpairs = n_part / 2
    density = 1.5736
    time_step = 0.001823
    temp = 298.0
    gamma = 20.0
    #l_bjerrum = 0.885^2 * e^2/(4*pi*epsilon_0*k_B*T)
    l_bjerrum = 130878.0 / temp

    num_steps_equilibration = 3000
    num_configs = 50
    integ_steps_per_config = 500

# Particle parameters
    types = {"Cl":          0, "Na": 1}
    numbers = {"Cl": n_ionpairs, "Na": n_ionpairs}
    charges = {"Cl": -1.0, "Na": 1.0}
    lj_sigmas = {"Cl":       3.85, "Na": 2.52}
    lj_epsilons = {"Cl":     192.45, "Na": 17.44}

    lj_cuts = {"Cl":  3.0 * lj_sigmas["Cl"],
               "Na": 3.0 * lj_sigmas["Na"]}

    masses = {"Cl":        35.453, "Na": 22.99}

# Setup System
    box_l = (n_ionpairs * sum(masses.values()) / density)**(1. / 3.)
    system.box_l = [box_l, box_l, box_l]
    system.periodicity = [1, 1, 1]
    system.time_step = time_step
    system.cell_system.skin = 0.3
    system.thermostat.set_langevin(kT=temp, gamma=gamma)

# Place particles
    q = 1
    l = box_l / n_ppside
    for i in range(n_ppside):
        for j in range(n_ppside):
            for k in range(n_ppside):
                p = numpy.array([i, j, k]) * l
                if q < 0:
                    system.part.add(id=len(
                        system.part), type=types["Cl"],  pos=p, q=charges["Cl"], mass=masses["Cl"])
                else:
                    system.part.add(id=len(
                        system.part), type=types["Na"],  pos=p, q=charges["Na"], mass=masses["Na"])

                q *= -1
            q *= -1
        q *= -1

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

# Lennard-Jones interactions parameters
    for s in [["Cl", "Na"], ["Cl", "Cl"], ["Na", "Na"]]:
        lj_sig = combination_rule_sigma(
            "Berthelot", lj_sigmas[s[0]], lj_sigmas[s[1]])
        lj_cut = combination_rule_sigma(
            "Berthelot", lj_cuts[s[0]], lj_cuts[s[1]])
        lj_eps = combination_rule_epsilon(
            "Lorentz", lj_epsilons[s[0]], lj_epsilons[s[1]])

        system.non_bonded_inter[types[s[0]], types[s[1]]].lennard_jones.set_params(
            epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")

    print("\n--->Tuning Electrostatics")
    #p3m = electrostatics.P3M(bjerrum_length=l_bjerrum, accuracy=1e-2, mesh=[84,84,84], cao=6)
    p3m = electrostatics.P3M(bjerrum_length=l_bjerrum, accuracy=1e-2)
    system.actors.add(p3m)

    print("\n--->Temperature Equilibration")
    system.time = 0.0
    for i in range(int(num_steps_equilibration / 100)):
        energy = system.analysis.energy()
        temp_measured = energy['kinetic'] / ((3.0 / 2.0) * n_part)
        print("t={0:.1f}, E_total={1:.2f}, E_coulomb={2:.2f}, T_cur={3:.4f}".format(system.time,
                                                                                    energy['total'],
                                                                                    energy['coulomb'],
                                                                                    temp_measured))
        system.integrator.run(100)
        visualizer.update()

    print("\n--->Integration")
    system.time = 0.0
    while (True):
        system.integrator.run(1)
        visualizer.update()


# Start simulation in seperate thread
t = Thread(target=main)
t.daemon = True
t.start()

# Start blocking visualizer
visualizer.start()
