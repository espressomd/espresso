include "myconfig.pxi"

from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "reaction_ensemble.hpp" namespace "ReactionEnsemble":

    ctypedef struct single_reaction:
        int* reactant_types
        int len_reactant_types
        int* reactant_coefficients
        int* product_types
        int len_product_types
        int* product_coefficients
        double equilibrium_constant
        int nu_bar


    ctypedef struct reaction_system:
        int nr_single_reactions
        single_reaction** reactions
        int* type_index
        int nr_different_types
        double* charges_of_types
        double standard_pressure_in_simulation_units
        double temperature_reaction_ensemble
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

    ctypedef struct collective_variable:
        double CV_minimum
        double CV_maximum
        double delta_CV
        double (*determine_current_state_in_collective_variable_with_index) (int, void*)
        int* corresponding_acid_types
        int nr_corresponding_acid_types
        int associated_type
        char* energy_boundaries_filename

    ctypedef struct wang_landau_system:
        int* histogram
        int len_histogram
        double* wang_landau_potential
        int nr_collective_variables
        collective_variable** collective_variables
        int* nr_subindices_of_collective_variable
        double wang_landau_parameter
        double initial_wang_landau_parameter
        int int_fill_value
        double double_fill_value
        int number_of_monte_carlo_moves_between_check_of_convergence
        double final_wang_landau_parameter
        int used_bins
        int monte_carlo_trial_moves
        int wang_landau_steps
        char* output_filename
        double* minimum_energies_at_flat_index
        double* maximum_energies_at_flat_index
        bool do_energy_reweighting
        int polymer_start_id
        int polymer_end_id
        bool fix_polymer
        bool do_not_sample_reaction_partition_function

    cdef cppclass c_reaction_ensemble "ReactionEnsemble::ReactionEnsemble":
        reaction_system m_current_reaction_system
        int do_reaction(int reaction_steps) except +
        bool do_global_mc_move_for_particles_of_type(int type, int start_id_polymer, int end_id_polymer, int particle_number_of_type, bool use_wang_landau)
        int find_index_of_type(int type) except +
        void set_cuboid_reaction_ensemble_volume()
        int check_reaction_ensemble() except +
        int m_accepted_configurational_MC_moves
        int m_tried_configurational_MC_moves
        int delete_particle (int p_id)
        void add_reaction(double equilibrium_constant, vector[int] _reactant_types, vector[int] _reactant_coefficients, vector[int] _product_types, vector[int] _product_coefficients) except +

    #///////////////////////////////////////////// Wang-Landau reaction ensemble algorithm
        wang_landau_system m_current_wang_landau_system
        void add_new_CV_degree_of_association(int associated_type, double CV_minimum, double CV_maximum, vector[int] _corresponding_acid_types)
        void add_new_CV_potential_energy(string filename, double delta_CV)
        int do_reaction_wang_landau() except +
        int update_maximum_and_minimum_energies_at_current_state()
        void write_out_preliminary_energy_run_results(char* filename)
        int write_wang_landau_checkpoint(char* identifier)
        int load_wang_landau_checkpoint(char* identifier)
        void write_wang_landau_results_to_file(char* full_path_to_output_filename)
        bool do_HMC_move_wang_landau()

    #///////////////////////////////////////////// Constant-pH ensemble
        double m_constant_pH
        void set_pH(double pH)
        int do_reaction_constant_pH() except +


cdef class reaction_ensemble:
    cdef _params
