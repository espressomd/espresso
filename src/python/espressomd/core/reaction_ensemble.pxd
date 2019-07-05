include "myconfig.pxi"

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp.string cimport string
from libcpp.map cimport map

cdef extern from "reaction_ensemble.hpp" namespace "ReactionEnsemble":

    ctypedef struct SingleReaction:
        vector[int] reactant_types
        vector[int] reactant_coefficients
        vector[int] product_types
        vector[int] product_coefficients
        double gamma
        int nu_bar
        double get_acceptance_rate()

    cdef cppclass CReactionAlgorithm "ReactionEnsemble::ReactionAlgorithm":
        int CReactionAlgorithm(int seed)
        int do_reaction(int reaction_steps) except +
        bool do_global_mc_move_for_particles_of_type(int type, int particle_number_of_type, bool use_wang_landau)
        void set_cuboid_reaction_ensemble_volume()
        int check_reaction_ensemble() except +
        double get_acceptance_rate_configurational_moves()
        int delete_particle(int p_id)
        void add_reaction(double gamma, vector[int] _reactant_types, vector[int] _reactant_coefficients, vector[int] _product_types, vector[int] _product_coefficients) except +

        vector[SingleReaction] reactions
        int nr_different_types
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

    cdef cppclass CReactionEnsemble "ReactionEnsemble::ReactionEnsemble"(CReactionAlgorithm):
        CReactionEnsemble(int seed)

    cdef cppclass CWangLandauReactionEnsemble "ReactionEnsemble::WangLandauReactionEnsemble"(CReactionAlgorithm):
        CWangLandauReactionEnsemble(int seed)
        double wang_landau_parameter
        double initial_wang_landau_parameter
        int number_of_monte_carlo_moves_between_check_of_convergence
        double final_wang_landau_parameter
        string output_filename
        vector[double] minimum_energies_at_flat_index
        vector[double] maximum_energies_at_flat_index
        bool do_not_sample_reaction_partition_function
        void add_new_CV_degree_of_association(int associated_type, double CV_minimum, double CV_maximum, vector[int] corresponding_acid_types)
        void add_new_CV_potential_energy(string filename, double delta_CV)
        int update_maximum_and_minimum_energies_at_current_state()
        void write_out_preliminary_energy_run_results(string filename)
        int write_wang_landau_checkpoint(string identifier)
        int load_wang_landau_checkpoint(string identifier)
        void write_wang_landau_results_to_file(string full_path_to_output_filename)

    cdef cppclass CConstantpHEnsemble "ReactionEnsemble::ConstantpHEnsemble"(CReactionAlgorithm):
        CConstantpHEnsemble(int seed)
        double m_constant_pH

    cdef cppclass CWidomInsertion "ReactionEnsemble::WidomInsertion"(CReactionAlgorithm):
        CWidomInsertion(int seed)
        pair[double, double] measure_excess_chemical_potential(int reaction_id)
