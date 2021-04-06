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
        bool do_global_mc_move_for_particles_of_type(int type, int particle_number_of_type, bool use_wang_landau)
        void set_cuboid_reaction_ensemble_volume()
        int check_reaction_method() except +
        double get_acceptance_rate_configurational_moves()
        int delete_particle(int p_id)
        void add_reaction(double gamma, vector[int] reactant_types, vector[int] reactant_coefficients, vector[int] product_types, vector[int] product_coefficients) except +
        void delete_reaction(int reaction_id)

        vector[SingleReaction] reactions
        map[int, double] charges_of_types
        double temperature
        double exclusion_radius
        double volume
        bool box_is_cylindric_around_z_axis
        double cyl_radius
        double cyl_x
        double cyl_y
        bool box_has_wall_constraints
        double slab_start_z
        double slab_end_z
        int non_interacting_type

cdef extern from "reaction_methods/ReactionEnsemble.hpp" namespace "ReactionMethods":

    cdef cppclass CReactionEnsemble "ReactionMethods::ReactionEnsemble"(CReactionAlgorithm):
        CReactionEnsemble(int seed)

cdef extern from "reaction_methods/WangLandauReactionEnsemble.hpp" namespace "ReactionMethods":

    cdef cppclass CWangLandauReactionEnsemble "ReactionMethods::WangLandauReactionEnsemble"(CReactionAlgorithm):
        CWangLandauReactionEnsemble(int seed)
        double wang_landau_parameter
        double final_wang_landau_parameter
        string output_filename
        vector[double] minimum_energies_at_flat_index
        vector[double] maximum_energies_at_flat_index
        bool do_not_sample_reaction_partition_function
        void add_new_CV_degree_of_association(int associated_type, double CV_minimum, double CV_maximum, vector[int] corresponding_acid_types)
        void add_new_CV_potential_energy(string filename, double delta_CV) except +
        int update_maximum_and_minimum_energies_at_current_state()
        void write_out_preliminary_energy_run_results(string filename) except +
        void write_wang_landau_checkpoint(string identifier) except +
        void load_wang_landau_checkpoint(string identifier) except +
        void write_wang_landau_results_to_file(string filename) except +

cdef extern from "reaction_methods/ConstantpHEnsemble.hpp" namespace "ReactionMethods":

    cdef cppclass CConstantpHEnsemble "ReactionMethods::ConstantpHEnsemble"(CReactionAlgorithm):
        CConstantpHEnsemble(int seed)
        double m_constant_pH

cdef extern from "reaction_methods/WidomInsertion.hpp" namespace "ReactionMethods":

    cdef cppclass CWidomInsertion "ReactionMethods::WidomInsertion"(CReactionAlgorithm):
        CWidomInsertion(int seed)
        pair[double, double] measure_excess_chemical_potential(int reaction_id) except +
