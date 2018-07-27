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
from espressomd import assert_features, electrostatics
import numpy

assert_features(["ELECTROSTATICS", "LENNARD_JONES"])

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
types = {"Anion":          0, "Cation": 1}
numbers = {"Anion": n_ionpairs, "Cation": n_ionpairs}
charges = {"Anion": -1.0, "Cation": 1.0}
lj_sigmas = {"Anion":        1.0, "Cation": 1.0}
lj_epsilons = {"Anion":        1.0, "Cation": 1.0}

WCA_cut = 2.**(1. / 6.)
lj_cuts = {"Anion":  WCA_cut * lj_sigmas["Anion"],
           "Cation": WCA_cut * lj_sigmas["Cation"]}

# Setup System
box_l = (n_part / density)**(1. / 3.)
system = espressomd.System(box_l=[box_l] * 3)
system.seed  = system.cell_system.get_state()['n_nodes'] * [1234]
system.periodicity = [1, 1, 1]
system.time_step = time_step
system.cell_system.skin = 0.3
system.thermostat.set_langevin(kT=temp, gamma=gamma)

# Place particles
for i in range(int(n_ionpairs)):
    system.part.add(id=len(system.part), type=types["Anion"],  pos=numpy.random.random(
        3) * box_l, q=charges["Anion"])
for i in range(int(n_ionpairs)):
    system.part.add(id=len(system.part), type=types["Cation"], pos=numpy.random.random(
        3) * box_l, q=charges["Cation"])


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
    lj_sig = combination_rule_sigma(
        "Berthelot", lj_sigmas[s[0]], lj_sigmas[s[1]])
    lj_cut = combination_rule_sigma("Berthelot", lj_cuts[s[0]], lj_cuts[s[1]])
    lj_eps = combination_rule_epsilon(
        "Lorentz", lj_epsilons[s[0]], lj_epsilons[s[1]])

    system.non_bonded_inter[types[s[0]], types[s[1]]].lennard_jones.set_params(
        epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift="auto")


print("\n--->Lennard Jones Equilibration")
max_sigma = max(lj_sigmas.values())
min_dist = 0.0
cap = 10.0
# Warmup Helper: Cold, highly damped system
system.thermostat.set_langevin(kT=temp * 0.1, gamma=gamma * 50.0)

while min_dist < max_sigma:
    # Warmup Helper: Cap max. force, increase slowly for overlapping particles
    min_dist = system.analysis.mindist([types["Anion"], types["Cation"]], [
                                       types["Anion"], types["Cation"]])
    cap += min_dist
# print min_dist, cap
    system.force_cap = cap
    system.integrator.run(10)

# Don't forget to reset thermostat, timestep and force cap
system.thermostat.set_langevin(kT=temp, gamma=gamma)
system.force_cap = 0

print("\n--->Tuning Electrostatics")
p3m = electrostatics.P3M(bjerrum_length=l_bjerrum, accuracy=1e-3)
system.actors.add(p3m)

print("\n--->Temperature Equilibration")
system.time = 0.0
for i in range(int(num_steps_equilibration / 100)):
    temp_measured = system.analysis.energy(
    )['kinetic'] / ((3.0 / 2.0) * n_part)
    print("t={0:.1f}, E_total={1:.2f}, E_coulomb={2:.2f}, T_cur={3:.4f}".format(system.time,
                                                                                system.analysis.energy()[
                                                                                    'total'],
                                                                                system.analysis.energy()[
                                                                                    'coulomb'],
                                                                                temp_measured))
    system.integrator.run(100)

print("\n--->Integration")
system.time = 0.0
temp_measured = []
for i in range(num_configs):
    temp_measured.append(system.analysis.energy()[
                         'kinetic'] / ((3.0 / 2.0) * n_part))
    print("t={0:.1f}, E_total={1:.2f}, E_coulomb={2:.2f}, T_cur={3:.4f}".format(system.time,
                                                                                system.analysis.energy()[
                                                                                    'total'],
                                                                                system.analysis.energy()[
                                                                                    'coulomb'],
                                                                                temp_measured[-1]))
    system.integrator.run(integ_steps_per_config)

    # Interally append particle configuration
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
