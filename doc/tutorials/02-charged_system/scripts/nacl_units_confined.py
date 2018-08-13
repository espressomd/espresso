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
from espressomd import electrostatics, electrostatic_extensions, assert_features
from espressomd.shapes import Wall
import numpy

assert_features(["ELECTROSTATICS", "CONSTRAINTS", "MASS", "LENNARD_JONES"])

system = espressomd.System(box_l=[1.0, 1.0, 1.0])
system.seed  = system.cell_system.get_state()['n_nodes'] * [1234]
numpy.random.seed(system.seed)

print("\n--->Setup system")

# System parameters
n_part = 500
n_ionpairs = n_part / 2
density = 1.1138
time_step = 0.001823
temp = 1198.3
gamma = 50
#l_bjerrum = 0.885^2 * e^2/(4*pi*epsilon_0*k_B*T)
l_bjerrum = 130878.0 / temp
#[E]=k_b*K/e/A
#Ez = U[V]/(8.61733e-5*box_l)
# Ez=465.116 -> 1V
Ez = 465.116

num_steps_equilibration = 3000
num_configs = 500
integ_steps_per_config = 100

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


energy = system.analysis.energy()
print("Before Minimization: E_total=", energy['total'])
system.minimize_energy.init(
    f_max=10, gamma=50.0, max_steps=1000, max_displacement=0.2)
system.minimize_energy.minimize()
energy = system.analysis.energy()
print("After Minimization: E_total=", energy['total'])

print("\n--->Tuning Electrostatics")
p3m = electrostatics.P3M(bjerrum_length=l_bjerrum, accuracy=1e-2)
system.actors.add(p3m)
elc = electrostatic_extensions.ELC(gap_size=elc_gap, maxPWerror=1e-3)
system.actors.add(elc)

for p in system.part:
    p.ext_force = [0, 0, Ez * p.q]

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


print("\n--->Integration")
bins = 100
z_dens_na = numpy.zeros(bins)
z_dens_cl = numpy.zeros(bins)
system.time = 0.0
cnt = 0

for i in range(num_configs):
    temp_measured = system.analysis.energy(
    )['kinetic'] / ((3.0 / 2.0) * n_part)
    print("t={0:.1f}, E_total={1:.2f}, E_coulomb={2:.2f}, T_cur={3:.4f}".format(system.time,
                                                                                system.analysis.energy()[
                                                                                    'total'],
                                                                                system.analysis.energy()[
                                                                                    'coulomb'],
                                                                                temp_measured))
    system.integrator.run(integ_steps_per_config)

    for p in system.part:
        bz = int(p.pos[2] / box_z * bins)
        if p.type == types["Na"]:
            z_dens_na[bz] += 1.0
        elif p.type == types["Cl"]:
            z_dens_cl[bz] += 1.0
    cnt += 1

print("\n--->Analysis")
# Average / Normalize with Volume
z_dens_na /= (cnt * box_volume / bins)
z_dens_cl /= (cnt * box_volume / bins)
z_values = numpy.linspace(0, box_l, num=bins)
res = numpy.column_stack((z_values, z_dens_na, z_dens_cl))
numpy.savetxt("z_density.data", res, header="#z rho_na(z) rho_cl(z)")
print("\n--->Written z_density.data")
print("\n--->Done")
