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
import espressomd
from espressomd import thermostat
from espressomd import electrostatics
from espressomd import code_info
from espressomd import integrate
import numpy
import cPickle as pickle

# Set system parameters
n_part = 200
density = 0.7
box_l = (n_part / density)**(1. / 3.)

# Select electrostatics method
method = "p3m"
#method = "memd"
implemented_methods = ["p3m"]
#implemented_methods = ["p3m", "memd"]

# Setup system geometry in Espresso
system = espressomd.System()

system.box_l = [box_l, box_l, box_l]
system.periodicity = [1, 1, 1]

# Place particles
q = 1.
ptype = 0

for i in xrange(n_part):
    #posx, posy, posz = box_l*numpy.random.random(3)
    q *= -1
    ptype = 1 - ptype
    system.part.add(id=i, type=ptype, pos=numpy.random.random(3)
                    * system.box_l, q=q)

# Simulation parameters
system.time_step = 0.01
system.skin = 0.3

# Thermostat
temp = 1.
gamma = 1.
thermostat.Thermostat().set_langevin(kT=temp, gamma=gamma)

# Lennard-Jones interactions
lj_sig = 1.
lj_eps = 1.
lj_cut = 2.**(1. / 6.)
system.non_bonded_inter[0, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")
system.non_bonded_inter[1, 0].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")
system.non_bonded_inter[1, 1].lennard_jones.set_params(
    epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")

# Check if electrostatics method is clear...
if method not in implemented_methods:
    print("Please select an electrostatics method in the script.")
    exit()

# Distinguish between different methods
if method == "p3m":
    p3m = electrostatics.P3M(bjerrum_length=10.0, accuracy=1e-3)
    system.actors.add(p3m)
    print("P3M parameter:\n")
    p3m_params = p3m.get_params()
    for key in p3m_params.keys():
        print("{} = {}".format(key, p3m_params[key]))
# elif method == "memd":
    # TODO

if "ROTATION" in code_info.features():
    deg_free = 6
else:
    deg_free = 3

# Integration parameters
integ_steps = 200
int_n_times = 20

# Warmup integration loop
print("\nWarmup")
for cap in xrange(20, 200, 20):
    print("t={0}, E={1}".format(system.time,
                                system.analysis.energy()['total']))
    system.non_bonded_inter.set_force_cap(cap)
    integrate.integrate(integ_steps)
system.non_bonded_inter.set_force_cap(0)

# Pickle system properties
with open("system_config", "w") as system_config:
    pickle.dump(system, system_config, -1)

# Main integration loop
print("\nEquilibration")
for i in xrange(int_n_times):
    temp = system.analysis.energy()['ideal'] / ((deg_free / 2.0) * n_part)
    print("t={0}, E={1}, T={2}".format(system.time,
                                       system.analysis.energy()['total'], temp))
    integrate.integrate(integ_steps)

    # Pickle particle data
    with open("config_{}".format(i), "w") as configfile:
        pickle.dump(system.part, configfile, -1)
