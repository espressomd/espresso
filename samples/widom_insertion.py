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
Measure the excess chemical potential of a charged WCA fluid via Widom's
insertion method.
"""
import numpy as np
import argparse

import espressomd
from espressomd import reaction_ensemble
from espressomd import electrostatics

required_features = ["WCA", "P3M"]
espressomd.assert_features(required_features)

parser = argparse.ArgumentParser(epilog=__doc__)
parser.add_argument('cs_bulk', type=float,
                    help="bulk salt concentration [1/sigma^3]")
args = parser.parse_args()

# System parameters
#############################################################
cs_bulk = args.cs_bulk
N0 = 70
box_l = (N0 / cs_bulk)**(1.0 / 3.0)

# Integration parameters
#############################################################
system = espressomd.System(box_l=[box_l, box_l, box_l])
np.random.seed(seed=42)
system.time_step = 0.01
system.cell_system.skin = 0.4
temperature = 1.0


#############################################################
#  Setup System                                             #
#############################################################

# Particle setup
#############################################################
# type 0 = HA
# type 1 = A-
# type 2 = H+

for i in range(N0):
    system.part.add(id=i, pos=np.random.random(3) * system.box_l, type=1, q=-1)
for i in range(N0, 2 * N0):
    system.part.add(id=i, pos=np.random.random(3) * system.box_l, type=2, q=1)

wca_eps = 1.0
wca_sig = 1.0
types = [0, 1, 2]
for type_1 in types:
    for type_2 in types:
        system.non_bonded_inter[type_1, type_2].wca.set_params(
            epsilon=wca_eps, sigma=wca_sig)

p3m = electrostatics.P3M(prefactor=2.0, accuracy=1e-3)
system.actors.add(p3m)
p3m_params = p3m.get_params()
for key, value in p3m_params.items():
    print("{} = {}".format(key, value))

# Warmup
#############################################################
# warmup integration (steepest descent)
warm_steps = 20
warm_n_times = 20
min_dist = 0.9 * wca_sig

# minimize energy using min_dist as the convergence criterion
system.integrator.set_steepest_descent(f_max=0, gamma=1e-3,
                                       max_displacement=0.01)
i = 0
while system.analysis.min_dist() < min_dist and i < warm_n_times:
    print("minimization: {:+.2e}".format(system.analysis.energy()["total"]))
    system.integrator.run(warm_steps)
    i += 1

print("minimization: {:+.2e}".format(system.analysis.energy()["total"]))
print()
system.integrator.set_vv()

# activate thermostat
system.thermostat.set_langevin(kT=temperature, gamma=1.0, seed=42)

widom = reaction_ensemble.WidomInsertion(
    temperature=temperature, seed=77)

# add insertion reaction
insertion_reaction_id = 0
widom.add_reaction(reactant_types=[],
                   reactant_coefficients=[], product_types=[1, 2],
                   product_coefficients=[1, 1], default_charges={1: -1, 2: +1})
print(widom.get_status())
system.setup_type_map([0, 1, 2])

n_iterations = 100
for i in range(n_iterations):
    for j in range(50):
        widom.measure_excess_chemical_potential(insertion_reaction_id)
    system.integrator.run(steps=500)
    if i % 20 == 0:
        print("mu_ex_pair ({:.4f}, +/- {:.4f})".format(
            *widom.measure_excess_chemical_potential(insertion_reaction_id)))
        print("HA", system.number_of_particles(type=0), "A-",
              system.number_of_particles(type=1), "H+",
              system.number_of_particles(type=2))

print("excess chemical potential for an ion pair ",
      widom.measure_excess_chemical_potential(insertion_reaction_id))
