#
# Copyright (C) 2013-2018 The ESPResSo project
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
from __future__ import division
import numpy as np

import espressomd
espressomd.assert_features(["ELECTROSTATICS", "LENNARD_JONES"])
from espressomd import code_info
from espressomd import analyze
from espressomd import integrate
from espressomd.interactions import *
from espressomd import reaction_ensemble
from espressomd import polymer
from espressomd.polymer import create_polymer
from espressomd import interactions
from espressomd import electrostatics
import sys

# System parameters
#
box_l = 56.134
l_bjerrum = 2.0
temperature = 1.

system = espressomd.System(box_l=[box_l] * 3)

# Integration parameters
#

system.set_random_state_PRNG()
np.random.seed(seed=system.seed)
system.time_step = 0.01
system.cell_system.skin = 10.  # only for tutorial purposes
system.cell_system.max_num_cells = 2744

system.thermostat.set_langevin(kT=temperature, gamma=1.0, seed=42)


#
# Setup System                                             #
#

# reaction method
mode = "reaction_ensemble"
# mode="constant_pH_ensemble"


# Particle setup
#
N_P = 1  # number of chains
MPC = 50  # monomers per chain
N0 = N_P * MPC  # total number of monomers
nNaOH = 0  # number of initial Na+OH-
nHCl = 0  # number of initial H+Cl- (additional H+'s)


type_HA = 0   # type 0 = HA
type_A = 1   # type 1 = A-
type_H = 2   # type 2 = H+
type_OH = 3   # type 3 = OH-
type_Na = 4   # type 4 = Na+
type_Cl = 5   # type 5 = Cl-

charges = {}
charges[type_HA] = 0
charges[type_A] = -1
charges[type_H] = 1
charges[type_OH] = -1
charges[type_Na] = 1
charges[type_Cl] = -1


# bonding interaction parameter
bond_l = 1.2  # bond length
kbond = 100  # force constant for harmonic bond
harmonic_bond = interactions.HarmonicBond(k=kbond, r_0=bond_l)
system.bonded_inter.add(harmonic_bond)

# non-bonding interactions (LJ)
lj_eps = 1.0
lj_sig = 1.0
lj_cut = 1.12246
lj_shift = 0.0


# setting up the polymer
polymer.create_polymer(
    N_P=N_P, bond_length=bond_l, MPC=MPC, start_id=0, bond=harmonic_bond,
    type_poly_neutral=type_HA, type_poly_charged=type_A, mode=0,
    val_poly=charges[type_A], start_pos=[0] * 3)
# setting up counterions
for i in range(N0):
    system.part.add(pos=np.random.random(3) *
                    system.box_l, type=type_H, q=charges[type_H])

# setting up other ions
# - Na+ and OH-
for i in range(nNaOH):
    system.part.add(pos=np.random.random(3) *
                    system.box_l, type=type_OH, q=charges[type_OH])
for i in range(nNaOH):
    system.part.add(pos=np.random.random(3) *
                    system.box_l, type=type_Na, q=charges[type_Na])
# - (additional) H+ and Cl-
for i in range(nHCl):
    system.part.add(pos=np.random.random(3) *
                    system.box_l, type=type_H, q=charges[type_H])
for i in range(nHCl):
    system.part.add(pos=np.random.random(3) *
                    system.box_l, type=type_Cl, q=charges[type_Cl])


# setting up LJ-interactions
for i in range(1, 5):
    for j in range(i + 1, 6):
        system.non_bonded_inter[i, j].lennard_jones.set_params(
            epsilon=lj_eps, sigma=lj_sig, cutoff=lj_cut, shift=lj_shift)


# Setting up electrostatics
#
p3m = electrostatics.P3M(prefactor=l_bjerrum * temperature, accuracy=1e-3)
system.actors.add(p3m)


K_diss = 0.002694
K_w = 10.**(-14) * 0.02694**2
RE = None
if(mode == "reaction_ensemble"):
    RE = reaction_ensemble.ReactionEnsemble(
        temperature=temperature, exclusion_radius=1)
elif(mode == "constant_pH_ensemble"):
    RE = reaction_ensemble.ConstantpHEnsemble(
        temperature=temperature, exclusion_radius=1)
    RE.constant_pH = 0

# HA <--> A- + H+
RE.add_reaction(
    gamma=K_diss, reactant_types=[type_HA], reactant_coefficients=[1],
    product_types=[type_A, type_H], product_coefficients=[1, 1],
    default_charges={type_HA: charges[type_HA], type_A: charges[type_A],
                     type_H: charges[type_H]})

# H2O autoprotolysis
RE.add_reaction(
    gamma=(1 / K_w), reactant_types=[type_H, type_OH],
    reactant_coefficients=[1, 1], product_types=[], product_coefficients=[],
    default_charges={type_H: charges[type_H], type_OH: charges[type_OH]})


print(RE.get_status())
system.setup_type_map([type_HA, type_A, type_H, type_OH, type_Na, type_Cl])


alpha = []
nHA = []
nA = []
nH = []
nOH = []
qdist = np.zeros(N0)
c = 0
n_iterations = 500  # this is for tutorial only, too few integration steps
n_steps_production = 10000
n_steps_thermalization = 2000

for i in range(n_steps_thermalization + n_steps_production):
    RE.reaction()
    system.integrator.run(n_iterations)
    print(i, ") HA", system.number_of_particles(type=type_HA), "A-",
          system.number_of_particles(type=type_A), "H+",
          system.number_of_particles(type=type_H), 'OH-',
          system.number_of_particles(type=type_OH), 'Cl-',
          system.number_of_particles(type=type_Cl), 'NA+',
          system.number_of_particles(type=type_Na))
    if (i > n_steps_thermalization):  # just a bit of thermalization before starting to gain informations about the properties of the system
        alpha.append(system.number_of_particles(type=type_A) / N0)
        nHA.append(system.number_of_particles(type=type_HA))
        nA.append(system.number_of_particles(type=type_A))
        nH.append(system.number_of_particles(type=type_H))
        nOH.append(system.number_of_particles(type=type_OH))

        c = c + 1

        for n in range(N0):
            qn = system.part[n].q
            qdist[n] = qdist[n] + qn
            print(qdist)


alpha_av = np.mean(alpha)
alpha_err = np.std(alpha) / np.sqrt(len(alpha))

nHA_av = np.mean(nHA)
nA_av = np.mean(nA)
nH_av = np.mean(nH)
nOH_av = np.mean(nOH)


print("\n<alpha> = {} (err = {})".format(alpha_av, alpha_err))
print("\n")
print("\n<nHA> = {} ".format(nHA_av))
print("\n<nA>  = {} ".format(nA_av))
print("\n<nH>  = {} ".format(nH_av))
print("\n<nOH> = {} ".format(nOH_av))


qdist = qdist / c
print('*******************************************')
print('*** charge distribution along the chain ***')
for i in range(N0):
    print(i + 1, qdist[i])
