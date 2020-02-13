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
Simulate two reacting monomers which are bonded via a harmonic potential.
The script aborts as soon as the abortion criterion in the Wang-Landau
algorithm is met: the Wang-Landau simulation runs until the Wang-Landau
potential has converged and then raises a Warning that it has converged,
effectively aborting the simulation.

With the setup of the Wang-Landau algorithm in this script, you sample the
density of states of a three-dimensional reacting harmonic oscillator as
a function of the two collective variables 1) degree of association and
2) potential energy.

The recorded Wang-Landau potential (which is updated during the simulation)
is written to the file :file:`WL_potential_out.dat`.

In this simulation setup the Wang-Landau potential is the density of states.
You can view the converged Wang-Landau potential e.g. via plotting with
gnuplot: ``splot "WL_potential_out.dat"``. As expected the three-dimensional
harmonic oscillator has a density of states which goes like
:math:`\\sqrt{E_{\\text{pot}}}`.

For a scientific description and different ways to use the algorithm please
consult https://pubs.acs.org/doi/full/10.1021/acs.jctc.6b00791
"""
import numpy as np

import espressomd
from espressomd import reaction_ensemble
from espressomd.interactions import HarmonicBond

# System parameters
#############################################################
box_l = 6 * np.sqrt(2)

# Integration parameters
#############################################################
system = espressomd.System(box_l=[box_l] * 3)
np.random.seed(seed=42)
system.time_step = 0.02
system.cell_system.skin = 0.4


#############################################################
#  Setup System                                             #
#############################################################


# Particle setup
#############################################################
# type 0 = HA
# type 1 = A-
# type 2 = H+

N0 = 1  # number of titratable units
K_diss = 0.0088

system.part.add(id=0, pos=[0, 0, 0] * system.box_l, type=3)
system.part.add(id=1, pos=[1.0, 1.0, 1.0] * system.box_l / 2.0, type=1)
system.part.add(id=2, pos=np.random.random() * system.box_l, type=2)
system.part.add(id=3, pos=np.random.random() * system.box_l, type=2)

# create a harmonic bond between the two reacting particles => the
# potential energy is quadratic in the elongation of the bond and
# therefore the density of states is known as the one of the harmonic
# oscillator
h = HarmonicBond(r_0=0, k=1)
system.bonded_inter[0] = h
system.part[0].add_bond((h, 1))


RE = reaction_ensemble.WangLandauReactionEnsemble(
    temperature=1, exclusion_radius=0, seed=77)
RE.add_reaction(gamma=K_diss, reactant_types=[0], reactant_coefficients=[1],
                product_types=[1, 2], product_coefficients=[1, 1],
                default_charges={0: 0, 1: -1, 2: +1})
print(RE.get_status())
system.setup_type_map([0, 1, 2, 3])


# initialize wang_landau
# generate preliminary_energy_run_results here, this should be done in a
# separate simulation without energy reweighting using the update energy
# functions
np.savetxt("energy_boundaries.dat", np.c_[[0, 1], [0, 0], [9, 9]],
           header="nbar E_min E_max")

RE.add_collective_variable_degree_of_association(
    associated_type=0, min=0, max=1, corresponding_acid_types=[0, 1])
RE.add_collective_variable_potential_energy(
    filename="energy_boundaries.dat", delta=0.05)
RE.set_wang_landau_parameters(
    final_wang_landau_parameter=1e-3,
    do_not_sample_reaction_partition_function=True,
    full_path_to_output_filename="WL_potential_out.dat")

i = 0
while True:
    RE.reaction()
    RE.displacement_mc_move_for_particles_of_type(3)
