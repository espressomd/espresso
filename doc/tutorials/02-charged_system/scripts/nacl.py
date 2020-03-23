#
# Copyright (C) 2010-2019 The ESPResSo project
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
from espressomd import assert_features, electrostatics
from espressomd.minimize_energy import steepest_descent
import numpy

assert_features(["ELECTROSTATICS", "WCA"])

print("\n--->Setup system")

# System parameters
n_part = 200
n_ionpairs = n_part / 2
density = 0.5
time_step = 0.01
temp = 1.0
gamma = 1.0
l_bjerrum = 7.0

num_steps_equilibration = 1000
num_configs = 100
integ_steps_per_config = 1000

# Particle parameters
types = {"Anion": 0, "Cation": 1}
numbers = {"Anion": n_ionpairs, "Cation": n_ionpairs}
charges = {"Anion": -1.0, "Cation": 1.0}
wca_sigmas = {"Anion": 1.0, "Cation": 1.0}
wca_epsilons = {"Anion": 1.0, "Cation": 1.0}

# Setup System
box_l = (n_part / density)**(1. / 3.)
system = espressomd.System(box_l=[box_l] * 3)
system.periodicity = [True, True, True]
system.time_step = time_step
system.cell_system.skin = 0.3

# Place particles
for i in range(int(n_ionpairs)):
    system.part.add(id=len(system.part), type=types["Anion"],
                    pos=numpy.random.random(3) * box_l, q=charges["Anion"])
for i in range(int(n_ionpairs)):
    system.part.add(id=len(system.part), type=types["Cation"],
                    pos=numpy.random.random(3) * box_l, q=charges["Cation"])


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
for s in [["Anion", "Cation"], ["Anion", "Anion"], ["Cation", "Cation"]]:
    wca_sig = combination_rule_sigma(
        "Berthelot", wca_sigmas[s[0]], wca_sigmas[s[1]])
    wca_eps = combination_rule_epsilon(
        "Lorentz", wca_epsilons[s[0]], wca_epsilons[s[1]])

    system.non_bonded_inter[types[s[0]], types[s[1]]].wca.set_params(
        epsilon=wca_eps, sigma=wca_sig)


print("\n--->WCA Equilibration")
max_sigma = max(wca_sigmas.values())
min_dist = 0.0

while min_dist < max_sigma:
    steepest_descent(system, f_max=0, gamma=10, max_steps=10,
                     max_displacement=max_sigma * 0.01)
    min_dist = system.analysis.min_dist()

# Set thermostat
system.thermostat.set_langevin(kT=temp, gamma=gamma, seed=42)

print("\n--->Tuning Electrostatics")
p3m = electrostatics.P3M(prefactor=l_bjerrum, accuracy=1e-3)
system.actors.add(p3m)

print("\n--->Temperature Equilibration")
system.time = 0.0
for i in range(int(num_steps_equilibration / 100)):
    temp_measured = system.analysis.energy()['kinetic'] / ((3. / 2.) * n_part)
    print("t={0:.1f}, E_total={1:.2f}, E_coulomb={2:.2f}, T_cur={3:.4f}"
          .format(system.time, system.analysis.energy()['total'],
                  system.analysis.energy()['coulomb'], temp_measured))
    system.integrator.run(100)

print("\n--->Integration")
system.time = 0.0
temp_measured = []
for i in range(num_configs):
    temp_measured.append(system.analysis.energy()['kinetic']
                         / ((3.0 / 2.0) * n_part))
    print("t={0:.1f}, E_total={1:.2f}, E_coulomb={2:.2f}, T_cur={3:.4f}"
          .format(system.time, system.analysis.energy()['total'],
                  system.analysis.energy()['coulomb'], temp_measured[-1]))
    system.integrator.run(integ_steps_per_config)

    # Internally append particle configuration
    system.analysis.append()


print("\n--->Analysis")
# Calculate the averaged rdfs
rdf_bins = 100
r_min = 0.0
r_max = system.box_l[0] / 2.0
r, rdf_00 = system.analysis.rdf(rdf_type='<rdf>',
                                type_list_a=[types["Anion"]],
                                type_list_b=[types["Anion"]],
                                r_min=r_min,
                                r_max=r_max,
                                r_bins=rdf_bins)

r, rdf_01 = system.analysis.rdf(rdf_type='<rdf>',
                                type_list_a=[types["Anion"]],
                                type_list_b=[types["Cation"]],
                                r_min=r_min,
                                r_max=r_max,
                                r_bins=rdf_bins)
# Write out the data
numpy.savetxt('rdf.data', numpy.c_[r, rdf_00, rdf_01])
print("\n--->Written rdf.data")
print("\n--->Done")
