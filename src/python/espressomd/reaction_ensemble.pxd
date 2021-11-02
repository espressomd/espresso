# Copyright (C) 2010-2019 The ESPResSo project
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
include "myconfig.pxi"

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.map cimport map
from .utils cimport Vector2d

cdef extern from "reaction_methods/SingleReaction.hpp" namespace "ReactionMethods":

    ctypedef struct SingleReaction:
        vector[int] reactant_types
        vector[int] reactant_coefficients
        vector[int] product_types
        vector[int] product_coefficients
        double gamma
        int nu_bar
        double get_acceptance_rate()

cdef extern from "reaction_methods/ReactionAlgorithm.hpp" namespace "ReactionMethods":

    cdef cppclass CReactionAlgorithm "ReactionMethods::ReactionAlgorithm":
        int CReactionAlgorithm(int seed)
        int do_reaction(int reaction_steps) except +
        bool do_global_mc_move_for_particles_of_type(int type, int particle_number_of_type)
        void set_cuboid_reaction_ensemble_volume()
        int check_reaction_method() except +
        double get_acceptance_rate_configurational_moves()
        int delete_particle(int p_id)
        void add_reaction(double gamma, vector[int] reactant_types, vector[int] reactant_coefficients, vector[int] product_types, vector[int] product_coefficients) except +
        void delete_reaction(int reaction_id)
        void set_cyl_constraint(double center_x, double center_y, double radius) except +
        void set_slab_constraint(double slab_start_z, double slab_end_z) except +
        void remove_constraint()
        Vector2d get_slab_constraint_parameters()

        vector[SingleReaction] reactions
        map[int, double] charges_of_types
        double kT
        double exclusion_radius
        double volume
        int non_interacting_type

cdef extern from "reaction_methods/ReactionEnsemble.hpp" namespace "ReactionMethods":

    cdef cppclass CReactionEnsemble "ReactionMethods::ReactionEnsemble"(CReactionAlgorithm):
        CReactionEnsemble(int seed)

cdef extern from "reaction_methods/ConstantpHEnsemble.hpp" namespace "ReactionMethods":

    cdef cppclass CConstantpHEnsemble "ReactionMethods::ConstantpHEnsemble"(CReactionAlgorithm):
        CConstantpHEnsemble(int seed)
        double m_constant_pH

cdef extern from "reaction_methods/WidomInsertion.hpp" namespace "ReactionMethods":

    cdef cppclass CWidomInsertion "ReactionMethods::WidomInsertion"(CReactionAlgorithm):
        CWidomInsertion(int seed)
        double calculate_particle_insertion_potential_energy(SingleReaction & current_reaction) except +
