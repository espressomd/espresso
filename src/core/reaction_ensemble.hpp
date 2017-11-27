#include "utils.hpp"
#include <string>
#include "energy.hpp"	//for calculate_current_potential_energy_of_system

namespace ReactionEnsemble {

struct single_reaction {
  // strict input to the algorithm
  std::vector<int> reactant_types;
  int len_reactant_types;
  std::vector<int> reactant_coefficients;
  std::vector<int> product_types;
  int len_product_types;
  std::vector<int> product_coefficients;
  double equilibrium_constant;
  // calculated values that are stored for performance reasons
  int nu_bar;
};

struct reaction_system {
  int nr_single_reactions;
  std::vector<single_reaction> reactions;
  std::vector<int> type_index;
  int nr_different_types; // is equal to length type_index
  std::vector<double> charges_of_types;
  double standard_pressure_in_simulation_units;
  double temperature;
  double exclusion_radius; // this is used as a kind of hard sphere radius, if
                           // particles are closer than that it is assumed that
                           // their interaction energy gets approximately
                           // infinite => these configurations do not contribute
                           // to the partition function and ensemble averages.
  double volume;
  bool box_is_cylindric_around_z_axis;
  double cyl_radius;
  double cyl_x;
  double cyl_y;
  bool box_has_wall_constraints;
  double slab_start_z;
  double slab_end_z;
  int non_interacting_type;
};

struct stored_particle_property {
  int p_id;
  double charge;
  int type;
};

struct collective_variable {
  double CV_minimum;
  double CV_maximum;
  double delta_CV;
  // for collective variables of type degree of association
  std::vector<int> corresponding_acid_types;
  int associated_type;
  // for collective variables of type energy
  std::string energy_boundaries_filename;
  
  double determine_current_state_in_current_collective_variable(){
  		if(!this->corresponding_acid_types.empty()){
			//found a collective variable which is not of the type of a degree_of_association association)
			return calculate_degree_of_association();
		}
		if(!this->energy_boundaries_filename.empty()){
			//found a collective variable which is not of the type of an energy
			return calculate_current_potential_energy_of_system();
		}else{
		    throw std::runtime_error("collective variable not implemented\n");
		}
  
  };
  
  private:
    /**
    * Returns the degree of association for the current collective variable. This is needed since you may use multiple degrees of association as collective variable for the Wang-Landau algorithm.
    */
    double calculate_degree_of_association(){
	    int total_number_of_corresponding_acid=0;
	    for(int corresponding_type_i=0; corresponding_type_i<this->corresponding_acid_types.size();corresponding_type_i++){
		    int num_of_current_type;
		    number_of_particles_with_type(this->corresponding_acid_types[corresponding_type_i],&num_of_current_type);
		    total_number_of_corresponding_acid+=num_of_current_type;
	    }
	    if(total_number_of_corresponding_acid==0){
	        throw std::runtime_error("Have you forgotten to specify all corresponding acid types? Total particle number of corresponding acid type is zero\n");
	    }
	    int num_of_associated_acid;
	    number_of_particles_with_type(this->associated_type,&num_of_associated_acid);
	    double degree_of_association=double(num_of_associated_acid)/total_number_of_corresponding_acid; //cast to double because otherwise any fractional part is lost
	    return degree_of_association;
    }
};

struct wang_landau_system {
  std::vector<int> histogram;
  std::vector<double> wang_landau_potential; // equals the logarithm to basis e of the
                                 // degeneracy of the states
  std::vector<collective_variable> collective_variables;
  std::vector<int> nr_subindices_of_collective_variable;
  double wang_landau_parameter; // equals the logarithm to basis e of the
                                // modification factor of the degeneracy of
                                // states when the state is visited
  double initial_wang_landau_parameter;

  int int_fill_value;
  double double_fill_value;

  int number_of_monte_carlo_moves_between_check_of_convergence;
  double final_wang_landau_parameter;
  int used_bins;               // for 1/t algorithm
  int monte_carlo_trial_moves; // for 1/t algorithm

  int wang_landau_steps; // may be used for performance improvements, when you
                         // do not want to record other observables in the tcl
                         // script
  std::string output_filename;

  std::vector<double>
      minimum_energies_at_flat_index; // only present in energy preparation run
  std::vector<double>
      maximum_energies_at_flat_index; // only present in energy preparation run

  bool do_energy_reweighting;
  int polymer_start_id;
  int polymer_end_id;
  bool fix_polymer;
  bool do_not_sample_reaction_partition_function;
};

class ReactionEnsemble {

public:
  ReactionEnsemble();
  ~ReactionEnsemble();

  reaction_system m_current_reaction_system = {
    0, // int nr_single_reactions
    std::vector<single_reaction>(), // reactions
    std::vector<int>(), //type_index
    0, // int nr_different_types
    std::vector<double>(), // charges_of_types
    -10.0, // double standard_pressure_in_simulation_units
    -10.0, // double temperature
    0.0, // double exclusion_radius
    -10, // double volume
    false, // bool box_is_cylindric_around_z_axis
    -10.0, // double cyl_radius
    -10.0, // double cyl_x
    -10.0, // double cyl_y
    false, // bool box_has_wall_constraints
    -10.0, // double slab_start_z
    -10.0, // double slab_end_z
    100 // int non_interacting_type
  };
  // the standard_pressure_in_simulation_units
  // is an input parameter for the reaction
  // ensemble;
    
  int m_accepted_configurational_MC_moves = 0;
  int m_tried_configurational_MC_moves = 0;
  bool m_system_is_in_1_over_t_regime = false;

  void set_cuboid_reaction_ensemble_volume();
  int do_reaction(int reaction_steps);
  void check_reaction_ensemble();
  int calculate_nu_bar(std::vector<int>& reactant_coefficients,
                       std::vector<int>& product_coefficients); // should only be used at when
                                               // defining a new reaction
  int update_type_index(std::vector<int>& reactant_types,
                        std::vector<int>& product_types); // assign different types an
                                                // index in a growing list that
                                                // starts at 0 and is
                                                // incremented by 1 for each new
                                                // type. the entry in the index
                                                // at place i is the
                                                // "type_value". therefore the
                                                // type of type "typevalue" has
                                                // the index i;
                                                bool generic_oneway_reaction(int reaction_id, int reaction_modus);
  int find_index_of_type(int type);
  bool do_global_mc_move_for_particles_of_type(int type, int start_id_polymer,
                                               int end_id_polymer,
                                               int particle_number_of_type,
                                               const bool use_wang_landau);
  int delete_particle(int p_id);
  void add_reaction(double equilibrium_constant,
                    std::vector<int> _reactant_types,
                    std::vector<int> _reactant_coefficients,
                    std::vector<int> _product_types,
                    std::vector<int> _product_coefficients);
  ///////////////////////////////////////////// Wang-Landau algorithm

  wang_landau_system m_current_wang_landau_system = {
    std::vector<int>(), // histogram
    std::vector<double>(), // wang_landau_potential
    std::vector<collective_variable>(), // collective_variables
    std::vector<int>(), // nr_subindices_of_collective_variable
    1.0, // double wang_landau_parameter
    1.0, // double initial_wang_landau_parameter
    -10, // int int_fill_value
    -10.0, // double double_fill_value
    5000, // int number_of_monte_carlo_moves_between_check_of_convergence
    0.00001, // double final_wang_landau_parameter
    -10, // int used_bins
    0, // int monte_carlo_trial_moves
    1, // int wang_landau_steps
    "", // output_filename
    std::vector<double>(), // minimum_energies_at_flat_index
    std::vector<double>(), // maximum_energies_at_flat_index
    false, // bool do_energy_reweighting
    -10, // int polymer_start_id
    -10, // int polymer_end_id
    false, // bool fix_polymer
    false // bool do_not_sample_reaction_partition_function
  };
  // use negative value as fill value since it cannot occur in
  // the wang_landau algorithm in the histogram and in the wang
  // landau potential, use only 1 wang_landau_steps if you want
  // to record other observables in the tcl script.

  void
  add_new_CV_degree_of_association(int associated_type, double CV_minimum,
                                   double CV_maximum,
                                   std::vector<int> _corresponding_acid_types);
  void add_new_CV_potential_energy(std::string filename, double delta_CV);
  int do_reaction_wang_landau();
  int update_maximum_and_minimum_energies_at_current_state(); // use for
                                                              // preliminary
                                                              // energy
                                                              // reweighting
                                                              // runs
  void write_out_preliminary_energy_run_results(std::string filename);

  double calculate_degree_of_association(int index_of_current_collective_variable);
  double calculate_current_potential_energy_of_system_wrap(int unimportant_int);

  // checkpointing, only designed to reassign values of a previous simulation to
  // a new simulation with the same initialization process
  int write_wang_landau_checkpoint(std::string identifier);
  int load_wang_landau_checkpoint(std::string identifier);
  void write_wang_landau_results_to_file(std::string full_path_to_output_filename);

  /////////////////////////////////////////////  Constant-pH Reactions
  int do_reaction_constant_pH();
  double m_constant_pH = -10;

private:
  // declaration of boring helper functions or member variables
  int m_invalid_charge =
      -10000; // this is the default charge which is assigned to a type which
              // occurs in a reaction. this charge has to be overwritten. if it
              // is not overwritten the reaction ensemble will complain.
  enum m_reaction_mode {
    reaction_ensemble_mode,
    reaction_ensemble_wang_landau_mode,
    constant_pH_mode
  };
  double factorial_Ni0_divided_by_factorial_Ni0_plus_nu_i(int Ni0, int nu_i);
  bool all_reactant_particles_exist(int reaction_id);
  int replace(int p_id, int desired_type);
  int create_particle(int desired_type);
  std::vector<int> m_empty_p_ids_smaller_than_max_seen_particle;
  int hide_particle(int p_id, int previous_type);
  void remove_bins_that_have_not_been_sampled();

    /**
    *Calculates the average of an integer array (used for the histogram of the Wang-Landau algorithm). It excludes values which are initialized to be negative. Those values indicate that the Wang-Landau algorithm should not sample those values. The values still occur in the list because we can only store "rectangular" value ranges.
    */

    template <typename T>
    double average_list_of_allowed_entries(T vector){
	    double result=0.0;
	    int counter_allowed_entries=0;
	    for(int i=0;i<vector.size();i++){
		    if(vector[i]>=0){ //checks for validity of index i (think of energy collective variables, in a cubic memory layout there will be indices which are not allowed by the energy boundaries. These values will be initalized with a negative fill value)
			    result+=(double) vector[i];
			    counter_allowed_entries+=1;
		    }
	    }
	    return result/counter_allowed_entries;
    }

    /**
    * Checks wether an integer is in an array of integers.
    */
    template <typename T>
    bool is_in_list(T value, std::vector<T> list){
	    for(int i=0;i<list.size();i++){
		    if(list[i]==value)
			    return true;
	    }
	    return false;
    }
    
  void append_particle_property_of_random_particle(
      int type, std::vector<stored_particle_property> &list_of_particles);
  void make_reaction_attempt(
      single_reaction& current_reaction,
      std::vector<stored_particle_property> &changed_particles_properties,
      std::vector<int> &p_ids_created_particles,
      std::vector<stored_particle_property> &hidden_particles_properties);
  double calculate_factorial_expression(single_reaction& current_reaction,
                                        int *old_particle_numbers);
  void restore_properties(std::vector<stored_particle_property>& property_list,
                          const int number_of_saved_properties);
  int add_types_to_index(std::vector<int>& type_list, int status_gc_init);
  double calculate_boltzmann_factor_reaction_ensemble(
      single_reaction& current_reaction, double E_pot_old, double E_pot_new,
      std::vector<int> &old_particle_numbers);
  void add_random_vector(double *vector, int len_vector,
                         double length_of_displacement);
  int get_flattened_index_wang_landau(
      std::vector<double>& current_state, std::vector<double>& collective_variables_minimum_values,
      std::vector<double>& collective_variables_maximum_values,
      std::vector<double>& delta_collective_variables_values, int nr_collective_variables);
  void get_random_position_in_box(double *out_pos);
  void
  get_random_position_in_box_enhanced_proposal_of_small_radii(double *out_pos);
  
  // declarations wang_landau
  int initialize_wang_landau(); // may first be called after all collective
                                // variables are added
  bool can_refine_wang_landau_one_over_t();
  bool achieved_desired_number_of_refinements_one_over_t();
  void refine_wang_landau_parameter_one_over_t();
  int get_flattened_index_wang_landau_of_current_state();
  void update_wang_landau_potential_and_histogram(
      int index_of_state_after_acceptance_or_rejection);
  double calculate_boltzmann_factor_reaction_ensemble_wang_landau(
      single_reaction& current_reaction, double E_pot_old, double E_pot_new,
      std::vector<int> &old_particle_numbers, int old_state_index,
      int new_state_index, bool only_make_configuration_changing_move);
  int get_flattened_index_wang_landau_without_energy_collective_variable(
      int flattened_index_with_energy_collective_variable,
      int collective_variable_index_energy_observable); // needed for energy
                                                        // collective variable
  double calculate_delta_degree_of_association(collective_variable& current_collective_variable);
  void initialize_histogram();
  void initialize_wang_landau_potential();
  void reset_histogram();
  int m_WL_accepted_moves = 0;
  int m_WL_tries = 0;

  // declarations constant pH
  int get_random_p_id();
  double
  calculate_boltzmann_factor_consant_pH(single_reaction& current_reaction,
                                        double E_pot_old, double E_pot_new);
};
}
